#!/usr/bin/make -f

PATH:=$(PATH):$(PROJECT_SOURCE_DIR):$(PROJECT_BINARY_DIR)/src/tool

TABLE=$(PROJECT_SOURCE_DIR)/brucella/$(TARGET).tsv

all: 11-$(TARGET)-resolve-blast.fasta 10-$(TARGET)-with-blast-pangenome2.fasta \
	06-$(TARGET)-pangenome1-rest.fasta 12-$(TARGET)-consensus.fasta

OP0=--debug --workers 2 --timing --seq-storage=compact
OP1=$(OP0) --export-alignment=1
OPN=--in-seqs 02-$(TARGET)-names.fasta
OP2=$(OP1) $(OPN)

01-$(TARGET)-orig.fasta: $(TABLE)
	get_seqs.py --table $(TABLE) --out $@

02-$(TARGET)-names.fasta: 01-$(TARGET)-orig.fasta
	replace_names.py --table $(TABLE) --fasta $< --out $@

03-$(TARGET)-anchors.fasta: 02-$(TARGET)-names.fasta
	find_anchors $(OP2) --out-file $@

04-$(TARGET)-pangenome1.fasta: 03-$(TARGET)-anchors.fasta
	make_pangenome $(OP2) --in-blocks $< --out-file $@

04a-$(TARGET)-pangenome1-stick.fasta: 04-$(TARGET)-pangenome1.fasta
	stick $(OP2) --in-blocks $< --out-file $@

05-$(TARGET)-pangenome1-with-rest.fasta: 04a-$(TARGET)-pangenome1-stick.fasta
	add_rest $(OP2) --in-blocks $< --out-file $@

06-$(TARGET)-pangenome1-rest.fasta: 04a-$(TARGET)-pangenome1-stick.fasta
	rest $(OP2) --in-blocks $< --out-file $@

07-$(TARGET)-aligned.fasta: 05-$(TARGET)-pangenome1-with-rest.fasta
	align_all $(OP2) --in-blocks $< --out-file $@

08-$(TARGET)-blast.fasta: 07-$(TARGET)-aligned.fasta
	blast_hits $(OP2) --in-blocks $< --out-file $@

09-$(TARGET)-with-blast.fasta: 07-$(TARGET)-aligned.fasta 08-$(TARGET)-blast.fasta
	cat $^ > $@

10-$(TARGET)-with-blast-pangenome2.fasta: 09-$(TARGET)-with-blast.fasta
	make_pangenome $(OP2) --in-blocks $< --out-file $@

11-$(TARGET)-resolve-blast.fasta: 07-$(TARGET)-aligned.fasta
	resolve_blast $(OP2) --in-blocks $< --out-file $@

12-$(TARGET)-consensus.fasta: 07-$(TARGET)-aligned.fasta
	consensus $(OP0) $(OPN) --in-blocks $< --out-consensus $@

