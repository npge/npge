#!/usr/bin/make -f

PATH=$${PATH}:$(PROJECT_SOURCE_DIR):$(PROJECT_BINARY_DIR)/src/tool

TABLE=$(PROJECT_SOURCE_DIR)/brucella/$(TARGET).tsv

all: 2-$(TARGET).fasta

1-$(TARGET).fasta:
	get_seqs.py --table $(TABLE) --out $@

2-$(TARGET).fasta: 1-$(TARGET).fasta
	replace_names.py --table $(TABLE) --fasta $< --out $@

