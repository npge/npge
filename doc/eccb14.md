# NPG-explorer: a new tool for nucleotide pangenome construction and analysis of closely related prokaryotic genomes

Boris Nagaev and Andrei Alexeevski

Genomes of closely related bacteria have highly similar
sequences of orthologous fragments but usually undergo
multiple rearrangements, long deletions, insertions of
mobile elements and occasionally horizontally transferred
regions.

We developed a new tool, Nucleotide PanGenome explorer
(NPG-explorer), designed for aligning and analysis of a
number of input closely related genomes. NPG-explorer
constructs nucleotide pangenome - a set of aligned blocks,
each block consisting of orthologous fragments. Fragments
having no orthologs are considered dummy blocks of one
fragment. Each nucleotide from input genomes belongs to
exactly one block of NPG (it is a reason for NPG
terminology). Minimum length of block (default 100 bp) and
minimum identity (default 90%) are algorithm parameters.
NPG-explorer iterates block detection algorithm until the
following criterion is satisfied: BLAST search
all-against-all block consensuses detects no hits of
appropriate size and identity.

In addition NPG-explorer provides: (1) Multiple alignments
of input chromosomes represented by a sequence of block
identifiers. These alignments allow to detect chromosomal
rearrangements. (2) File with consensus sequences of all
blocks and file with description of all mutations with
respect to consensuses. Thus, all input genome sequences
can be completely reconstructed from these two files. (3)
Phylogenetic trees of core blocks and of whole genomes.
Core blocks are those that contain exactly one fragment of
each genome. These trees are computed on the base of
diagnostic positions in block alignments. (4) All gene
annotations, mapped on blocks. This data are useful for
detection and correction mis-annotations, gene corruption
etc.

NPG-visualization tool presents interactively a list of
blocks, the alignment with mapped genes, alignments of
block identifiers.

NPG-explorer is written in C++ and is licensed under the
GNU GPL. Simple script language for program modules
invocation is introduced.

NPG-explorer was applied to 17 genomes of Brucella genus,
each genome consists of about 3 Mb. NPG-explorer worked
approximately one hour on Intel Core i5 computer. Among 527
detected blocks with two or more fragments there are 270
core blocks, they cover 95% of nucleotides. Average
sequence similarity in core blocks is 99.2%. Phylogenetic
tree of genomes computed by NPG-explorer by using
diagnostic positions is in agreement with published data
for 10 Brucella genomes [1]. Detected 25787 point mutations
and 2334 deletions of more than two bp describes the
evolution of sequences within blocks. The program found
large translocation from first to second chromosome in
Brucella suis ATCC 23445 and large inversion in chromosome
2 of Brucella abortus, also described earlier [2].

The work was supported by RFBR grants 14-04-01693,
13-07-00969.

    [1] Wattam et al., J.Bacteriology, 191:3569-79 (2009)
    [2] Tsoktouridis et al., J.Bacteriology, 185:6130-6 (2003)
