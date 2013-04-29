#!/usr/bin/make -f

PATH:=$(PATH):$(PROJECT_SOURCE_DIR):$(PROJECT_BINARY_DIR)/src/tool

TABLE=$(PROJECT_SOURCE_DIR)/brucella/$(TARGET).tsv

OP0=--debug --workers 2 --timing --seq-storage=compact
OP1=$(OP0) --export-alignment=1

01-$(TARGET)-orig.fasta: $(TABLE)
	get_seqs.py --table $(TABLE) --out $@

02-$(TARGET)-names.fasta: 01-$(TARGET)-orig.fasta
	replace_names.py --table $(TABLE) --fasta $< --out $@

03-$(TARGET)-anchors.fasta: 02-$(TARGET)-names.fasta
	find_anchors $(OP1) --in-blocks=$< --out-file $@

04-$(TARGET)-pangenome1.fasta: 03-$(TARGET)-anchors.fasta
	make_pangenome $(OP1) --in-blocks $< --out-file $@ --max-spreading 2.0

04a-$(TARGET)-pangenome1-stick.fasta: 04-$(TARGET)-pangenome1.fasta
	stick $(OP1) --in-blocks $< --out-file $@

05-$(TARGET)-pangenome1-with-rest.fasta: 04a-$(TARGET)-pangenome1-stick.fasta
	add_rest $(OP1) --in-blocks $< --out-file $@

06-$(TARGET)-pangenome1-rest.fasta: 04a-$(TARGET)-pangenome1-stick.fasta
	rest $(OP1) --in-blocks $< --out-file $@

07-$(TARGET)-aligned.fasta: 05-$(TARGET)-pangenome1-with-rest.fasta
	align_all $(OP0) --in-blocks $< --out-file $@

08-$(TARGET)-blast.fasta: 07-$(TARGET)-aligned.fasta
	blast_hits $(OP1) --in-blocks $< --out-file $@

09-$(TARGET)-with-blast.fasta: 07-$(TARGET)-aligned.fasta 08-$(TARGET)-blast.fasta
	cat $^ > $@

10-$(TARGET)-with-blast-cleanup.fasta: 09-$(TARGET)-with-blast.fasta
	make_pangenome $(OP1) --in-blocks $< --out-file $@

10a-$(TARGET)-pangenome2.fasta: 10-$(TARGET)-with-blast-cleanup.fasta
	add_rest $(OP1) --in-blocks $< --out-file $@

10b-$(TARGET)-pangenome2-aligned.fasta: 10a-$(TARGET)-pangenome2.fasta
	align_all $(OP0) --in-blocks $< --out-file $@

11-$(TARGET)-resolve-blast.fasta: 07-$(TARGET)-aligned.fasta
	resolve_blast $(OP1) --in-blocks $< --out-file $@

11a-$(TARGET)-resolve-blast-aligned.fasta: 11-$(TARGET)-resolve-blast.fasta
	align_all $(OP0) --in-blocks $< --out-file $@

12-$(TARGET)-consensus.fasta: 02-$(TARGET)-names.fasta 07-$(TARGET)-aligned.fasta
	consensus $(OP0) --in-blocks $^ --out-consensus $@

13-$(TARGET)-good-pangenome.fasta: 07-$(TARGET)-aligned.fasta
	blast_small_blocks.br $(OP0) --import-alignment=1 \
		--in-blocks $^ --out-file $@-temp1
	blast_small_blocks.br $(OP0) --import-alignment=1 \
		--in-blocks $@-temp1 --out-file $@-temp2
	gaps.br $(OP0) --import-alignment=1 --cut-gaps-mode=strict \
		--in-blocks $@-temp2 --out-file $@

13a-$(TARGET)-good-pangenome.fasta.ip: 13-$(TARGET)-good-pangenome.fasta
	meta_processor.br --processor IsPangenome $(OP0) --import-alignment=1 \
		--in-blocks $^ --out-file /dev/null > $@

13b-$(TARGET)-good-pangenome.fasta.stats: 13-$(TARGET)-good-pangenome.fasta
	stats $(OP0) --import-alignment=1 --in-blocks $^ > $@

13c-$(TARGET)-good-pangenome.fasta.bi: 13-$(TARGET)-good-pangenome.fasta
	blockinfo $(OP0) --import-alignment=1 --in-blocks $^ > $@

14-$(TARGET)-genes.txt:
	get_seqs.py $(OP0) --table $(TABLE) --out $@ --type genes

15-$(TARGET)-genes.fasta: 14-$(TARGET)-genes.txt 02-$(TARGET)-names.fasta
	add_genes.br --in-blocks 02-$(TARGET)-names.fasta \
		--in-genes 14-$(TARGET)-genes.txt --out-file $@

16-$(TARGET)-partition.fasta: 15-$(TARGET)-genes.fasta 13-$(TARGET)-good-pangenome.fasta
	partition.br $(OP0) --import-alignment=1 \
		--pangenome-in-blocks 13-$(TARGET)-good-pangenome.fasta \
		--genes-in-blocks 15-$(TARGET)-genes.fasta --out-file $@

16a-$(TARGET)-partition.xls: 15-$(TARGET)-genes.fasta 13-$(TARGET)-good-pangenome.fasta
	print_partition.br $(OP0) --import-alignment=1 \
		--pangenome-in-blocks 13-$(TARGET)-good-pangenome.fasta \
		--genes-in-blocks 15-$(TARGET)-genes.fasta --file $@

