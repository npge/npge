#!/usr/bin/make -f

PATH=$${PATH}:$(PROJECT_SOURCE_DIR):$(PROJECT_BINARY_DIR)/src/tool

TABLE=$(PROJECT_SOURCE_DIR)/brucella/$(TARGET).tsv

all: 4-$(TARGET)-pangenome1.fasta

OP1=--debug --workers 2 --timing --seq-storage=compact
OP2=$(OP1) --in-seqs 2-$(TARGET)-names.fasta

1-$(TARGET)-orig.fasta:
	get_seqs.py --table $(TABLE) --out $@

2-$(TARGET)-names.fasta: 1-$(TARGET)-orig.fasta
	replace_names.py --table $(TABLE) --fasta $< --out $@

3-$(TARGET)-anchors.fasta: 2-$(TARGET)-names.fasta
	find_anchors $(OP2) --out-file $@

4-$(TARGET)-pangenome1.fasta: 3-$(TARGET)-anchors.fasta
	make_pangenome $(OP2) --in-blocks $< --out-file $@

