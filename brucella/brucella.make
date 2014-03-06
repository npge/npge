#!/usr/bin/make -f

# bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
# Copyright (C) 2012 Boris Nagaev
#
# See the LICENSE file for terms of use.

PATH:=$(PATH):$(PROJECT_SOURCE_DIR):$(PROJECT_BINARY_DIR)/src/tool

TABLE=$(PROJECT_SOURCE_DIR)/brucella/$(TARGET).tsv

OP0=--debug --workers 2 --timing=1
OP1=$(OP0)

01-$(TARGET)-orig.fasta: $(TABLE)
	get_seqs.py --table $(TABLE) --out $@

02-$(TARGET)-names.fasta: 01-$(TARGET)-orig.fasta
	replace_names.py --table $(TABLE) --fasta $< --out $@

03-$(TARGET)-anchors.fasta: 02-$(TARGET)-names.fasta
	find_anchors.br $(OP1) --in-blocks=$< --out-file $@

13-$(TARGET)-good-pangenome.fasta: 02-$(TARGET)-names.fasta
	make_pangenome.br $(OP1) --in-blocks $< --out-file $@

13a-$(TARGET)-good-pangenome.fasta.ip: 13-$(TARGET)-good-pangenome.fasta
	meta_processor.br --processor IsPangenome $(OP0) \
		--in-blocks $^ --out-file /dev/null > $@

13b-$(TARGET)-good-pangenome.fasta.stats: 13-$(TARGET)-good-pangenome.fasta
	stats.br $(OP0) --in-blocks $^ > $@

13c-$(TARGET)-good-pangenome.fasta.bi: 13-$(TARGET)-good-pangenome.fasta
	blockinfo.br $(OP0) --in-blocks $^ > $@

14-$(TARGET)-genes.txt:
	get_seqs.py --table $(TABLE) --out $@ --type genes

15-$(TARGET)-genes.fasta: 14-$(TARGET)-genes.txt 02-$(TARGET)-names.fasta
	add_genes.br --in-blocks 02-$(TARGET)-names.fasta \
		--in-genes 14-$(TARGET)-genes.txt --out-file $@

16-$(TARGET)-partition.fasta: 15-$(TARGET)-genes.fasta 13-$(TARGET)-good-pangenome.fasta
	partition.br $(OP0) \
		--pangenome-in-blocks 13-$(TARGET)-good-pangenome.fasta \
		--genes-in-blocks 15-$(TARGET)-genes.fasta --out-file $@

16a-$(TARGET)-partition.xls: 15-$(TARGET)-genes.fasta 13-$(TARGET)-good-pangenome.fasta
	print_partition.br $(OP0) \
		--pangenome-in-blocks 13-$(TARGET)-good-pangenome.fasta \
		--genes-in-blocks 15-$(TARGET)-genes.fasta --file $@

17-$(TARGET)-gene-groups.fasta: 16-$(TARGET)-partition.fasta 13-$(TARGET)-good-pangenome.fasta
	find_gene_groups.br $(OP0) \
		--pangenome-in-blocks 13-$(TARGET)-good-pangenome.fasta \
		--genes-in-blocks 16-$(TARGET)-partition.fasta --out-file $@

17a-$(TARGET)-gene-parts.xls: 13-$(TARGET)-good-pangenome.fasta \
		16-$(TARGET)-partition.fasta 17-$(TARGET)-gene-groups.fasta
	print_gene_parts.br $(OP0) \
		--pangenome-in-blocks 13-$(TARGET)-good-pangenome.fasta \
		--genes-in-blocks 16-$(TARGET)-partition.fasta \
		--groups-in-blocks 17-$(TARGET)-gene-groups.fasta \
		--file $@

17b-$(TARGET)-gene-groups.xls: 13-$(TARGET)-good-pangenome.fasta \
		16-$(TARGET)-partition.fasta 17-$(TARGET)-gene-groups.fasta
	print_gene_groups.br $(OP0) \
		--pangenome-in-blocks 13-$(TARGET)-good-pangenome.fasta \
		--genes-in-blocks 16-$(TARGET)-partition.fasta \
		--groups-in-blocks 17-$(TARGET)-gene-groups.fasta \
		--file $@

20-$(TARGET)-good-pangenome-stem.fasta: 13-$(TARGET)-good-pangenome.fasta
	meta_processor.br --processor Stem $(OP0) \
		--in-blocks $^ --out-file $@ --skip-rest=1

20a-$(TARGET)-good-pangenome-stem.fasta.stats: 20-$(TARGET)-good-pangenome-stem.fasta
	stats.br $(OP0) --in-blocks $^ > $@

20b-$(TARGET)-good-pangenome-stem.fasta.bi: 20-$(TARGET)-good-pangenome-stem.fasta
	blockinfo.br $(OP0) --in-blocks $^ > $@

