#!/bin/bash

TARGET=$1
PROJECT_SOURCE_DIR=$2
PROJECT_BINARY_DIR=$3

export PATH=$PATH:$PROJECT_SOURCE_DIR:$PROJECT_BINARY_DIR/src/tool

TABLE=$PROJECT_SOURCE_DIR/brucella/$TARGET.tsv

get_seqs.py --table $TABLE --out 1-$TARGET.fasta
replace_names.py --table $TABLE --fasta 1-$TARGET.fasta --out 2-$TARGET.fasta

echo "Ok"
