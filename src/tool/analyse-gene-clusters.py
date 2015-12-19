#!/usr/bin/python

""" Analyze clusters of genes in output of find-genes-clusters.py

Reads file from stdin.
"""

import sys

clusters = {}

# read clusters file

for line in sys.stdin:
    line = line.strip()
    fields = line.split('\t')
    if fields[0] == "cluster_id":
        continue  # header
    cluster_id = int(fields[0])
    gene = fields[1]
    if cluster_id not in clusters:
        clusters[cluster_id] = {}
    clusters[cluster_id][gene] = []
    START_FIELD = 2
    GROUP_LENGTH = 5
    for start_field in range(START_FIELD, len(fields), GROUP_LENGTH):
        (block_name, block_min, block_max, block_ori, gene_min) = \
            tuple(fields[start_field : start_field + GROUP_LENGTH])
        if block_name:
            block_min = int(block_min)
            block_max = int(block_max)
            block_ori = int(block_ori)
            gene_min = int(gene_min)
        clusters[cluster_id][gene].append({
            'block_name': block_name,
            'block_min': block_min,
            'block_max': block_max,
            'block_ori': block_ori,
            'gene_min': gene_min,
        })

# find clusters with identical annotation

def endInBlock(part, end):
    ori = (end == 'start') and 1 or -1
    if part['block_ori'] * ori == 1:
        return part['block_min']
    else:
        return part['block_max']

def classify(genes):
    good = True
    same_stop = True
    first_parts = genes.values()[0]
    first_start_pos = endInBlock(first_parts[0], 'start')
    first_start_block = first_parts[0]['block_name']
    first_stop_pos = endInBlock(first_parts[-1], 'stop')
    first_stop_block = first_parts[-1]['block_name']
    first_tuple = (
        first_start_pos,
        first_start_block,
        first_stop_pos,
        first_stop_block,
    )
    first_tuple_stop = (
        first_stop_pos,
        first_stop_block,
    )
    for parts in genes.values():
        start_pos = endInBlock(parts[0], 'start')
        start_block = parts[0]['block_name']
        stop_pos = endInBlock(parts[-1], 'stop')
        stop_block = parts[-1]['block_name']
        this_tuple = (
            start_pos,
            start_block,
            stop_pos,
            stop_block,
        )
        this_tuple_stop = (
            stop_pos,
            stop_block,
        )
        if this_tuple != first_tuple:
            good = False
        if this_tuple_stop != first_tuple_stop:
            same_stop = False
    if good:
        return 'good'
    elif same_stop:
        return 'same_stop'
    else:
        return 'bad'

for cluster_id, genes in clusters.items():
    print("%d\t%s" % (cluster_id, classify(genes)))
