#!/usr/bin/python

""" Find clusters of genes in partition-grouped.tsv

Reads file partition-grouped.tsv from stdin, prints
output table "cluster_id", "fragment_id".

How to use:
$ find-genes-clusters.py 71 < genes/partition-grouped.tsv
"""

import sys
import logging

overlap_threshold = int(sys.argv[1])  # 71

def keepBlock(npg_block):
    return True

gene2blocks = {}

for line in sys.stdin:
    line = line.strip()
    fields = line.split('\t')
    sequence = fields[0]
    if sequence == 'sequence':
        continue  # header
    sequence_start = int(fields[1])
    sequence_stop = int(fields[2])
    gene = "%s_%d_%d" % (sequence, sequence_start, sequence_stop)
    assert gene not in gene2blocks
    gene2blocks[gene] = {}
    for block_field_index in range(4, len(fields), 6):
        npg_block = fields[block_field_index]
        npg_block_min = int(fields[block_field_index + 1])
        npg_block_max = int(fields[block_field_index + 2])
        npg_block_ori = int(fields[block_field_index + 3])
        gene_start = int(fields[block_field_index + 4])
        gene_stop = int(fields[block_field_index + 5])
        if not keepBlock(npg_block):
            continue
        if npg_block in gene2blocks[gene]:
            logging.warning("Gene %s overlaps block %s multiple times",
                    gene, npg_block)
        gene2blocks[gene][npg_block] = {
            "npg_block_min": npg_block_min,
            "npg_block_max": npg_block_max,
            "npg_block_ori": npg_block_ori,
            "gene_start": gene_start,
            "gene_stop": gene_stop,
        }
    if not gene2blocks[gene]:
        del gene2blocks[gene]

block2genes = {}

for gene, blocks in gene2blocks.items():
    for npg_block in blocks.keys():
        if npg_block not in block2genes:
            block2genes[npg_block] = set()
        block2genes[npg_block].add(gene)

genes2overlap = {}

for npg_block, genes in block2genes.items():
    for gene1 in genes:
        for gene2 in genes:
            if gene1 < gene2:
                block_min1 = gene2blocks[gene1][npg_block]["npg_block_min"]
                block_max1 = gene2blocks[gene1][npg_block]["npg_block_max"]
                block_min2 = gene2blocks[gene2][npg_block]["npg_block_min"]
                block_max2 = gene2blocks[gene2][npg_block]["npg_block_max"]
                max_min = max(block_min1, block_min2)
                min_max = min(block_max1, block_max2)
                if max_min <= min_max:
                    common = min_max - max_min + 1
                    fs = "%s\t%s" % (gene1, gene2)
                    if fs not in genes2overlap:
                        genes2overlap[fs] = 0
                    genes2overlap[fs] += common

graph = {}

for fs, common in genes2overlap.items():
    if common >= overlap_threshold:
        gene1, gene2 = list(fs.split())
        if gene1 not in graph:
            graph[gene1] = []
        graph[gene1].append(gene2)
        if gene2 not in graph:
            graph[gene2] = []
        graph[gene2].append(gene1)

cluster_id = 0
seen = set()

def findConnected(start_gene):
    connected_genes = []
    queue = set([start_gene])
    while queue:
        gene1 = queue.pop()
        if gene1 not in seen:
            seen.add(gene1)
            connected_genes.append(gene1)
            if gene1 in graph:
                for gene2 in graph[gene1]:
                    queue.add(gene2)
    return connected_genes

for start_gene in gene2blocks:
    if start_gene not in seen:
        cluster_id += 1
        cluster = list(findConnected(start_gene))
        cluster_blocks = set()
        for gene in cluster:
            cluster_blocks |= set(gene2blocks[gene].keys())
        def key(block1):
            for gene in cluster:
                block2props = gene2blocks[gene]
                if block1 in block2props:
                    return block2props[block1]["gene_start"]
        cluster_blocks = sorted(cluster_blocks, key=key)
        for gene in cluster:
            blocks = set(gene2blocks[gene].keys())
            in_blocks = []
            for block in cluster_blocks:
                # align blocks in cluster
                if block in blocks:
                    npg_block = block
                    npg_block_min = str(gene2blocks[gene][block]["npg_block_min"])
                    npg_block_max = str(gene2blocks[gene][block]["npg_block_max"])
                    npg_block_ori = str(gene2blocks[gene][block]["npg_block_ori"])
                    gene_start = str(gene2blocks[gene][block]["gene_start"])
                else:
                    npg_block = ''
                    npg_block_min = ''
                    npg_block_max = ''
                    npg_block_ori = ''
                    gene_start = ''
                in_blocks += [npg_block, npg_block_min, npg_block_max, npg_block_ori, gene_start]
            print("%d\t%s\t%s" %
                    (cluster_id, gene, '\t'.join(in_blocks)))
