#!/usr/bin/make -f

PATH=$${PATH}:$(PROJECT_SOURCE_DIR):$(PROJECT_BINARY_DIR)/src/tool

TABLE=$(PROJECT_SOURCE_DIR)/brucella/$(TARGET).tsv

all: 3-$(TARGET).fasta

OP1=--debug --workers 2 --timing --seq-storage=compact
OP2=$(OP1) --in-seqs 2-$(TARGET).fasta

1-$(TARGET).fasta:
	get_seqs.py --table $(TABLE) --out $@

2-$(TARGET).fasta: 1-$(TARGET).fasta
	replace_names.py --table $(TABLE) --fasta $< --out $@

3-$(TARGET).fasta: 2-$(TARGET).fasta
	find_anchors $(OP2) --out-file 3-$(TARGET).fasta

