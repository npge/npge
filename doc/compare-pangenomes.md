# Comparison of pangenomes

## Product of pangenomes

There are pangenomes A and B, built on the same set of sequences.

**Definition**: given pangenome C such that:

  1. every block from a pangenome C overlaps exactly one block
     from A and exactly one block from B and
  2. number of blocks in pangenome C is minimum

the pangenome C is called a **product** of pangenomes A and B.

### Algorithm of pangenome multiplication

**Goal**: build pangenome C which is a product of pangenomes A and B.

**Steps**. Let us describe the algorithm with pseudo-code:

```
FOR EACH block a FROM pangenome A:
  CREATE a dictionary D, mapping block from B to block from C;
  FOR EACH fragment f_a FROM block a:
    FIND fragment list ff_b from pangenome B, overlapping fragment f_a;
    FOR EACH fragment f_b FROM fragment list ff_b:
      FIND block b from pangenome B, containing f_b;
      IF dictionary D doesn't contain b:
        ADD pair (b, empty block) to dictionary D;
      FIND fragment f_c which is an intersection of f_a and f_b;
      ADD fragment f_c to block D[b];
  ADD all blocks from dictionary D to answer.
```

This algorithm can parallelized by block a.

### Properties of product of pangenomes

  * if pangenomes A and B are equal, then C=A=B;
  * if a block from A is equal to a block from B B, then a product
    of A and B contains the equal block;
  * if a block from A is similar to a block from B B, then a product
    of A and B contains the intersection of two blocks.

## Classification of blocks of product of pangenomes

There are pangenomes A and B built on the same set of sequences
and product C of pangenomes A and B.

Separate blocks of the product C to two types:

  * **common blocks**,
  * **conflict blocks**.

A block from pangenomes A and B can overlap with not more than one
common block. Other blocks from C are conflict blocks.
While choosing common blocks, blocks from product are ranked with
normal rules (by the number of fragments then by length).
Thus, better blocks tend to be marked as common blocks
and worse blocks (minor, unique) - as conflict blocks.

### Algorithm of classification of blocks of product

The algorithm of classification follows from the definition.

Make a set `used` which stores "used" blocks from A and B.
Sort blocks of the product with normal rules (by the number of
fragments then by length).
Iterate blocks from the best to the worst.
For each from c find corresponding block a in pangenome A and
block b in pangenome B.
Check if set `used` contains blocks a and b.
If it doesn't contain any of them, add blocks a and b to
set `used` and mark block c as a common block.
Otherwise mark block c as a conflict block.

### Properties of classification of blocks of product of pangenomes

  * if pangenomes A and B are equal, then all blocks from the product
    are common blocks;
  * if a block from A is equal to a block from B, then the equal block
    from the product is a common block;
  * is a block from A is similar to a block from B, then the intersection
    of these two blocks is a common block and other parts of them are
    conflict blocks.

## Distance between pangenomes

There are pangenomes A and B built on the same set of sequences,
the product C of pangenomes A and B and the classification of
the blocks of the product.

The **absolute distance** between pangenomes is the total number of
nucleotides in the conflict blocks excluding blocks corresponding
to minor blocks in both pangenomes.

The **relative distance** between pangenomes is the absolute distance
between the pangenomes divided by the sum of lengths of sequences
in a pangenome.

### Properties of distance between pangenomes

  * relative distance has a value from the interval [0, 1);
  * if pangenomes A and B are equal, then the distance is 0;
  * if a block from A is equal to a block from B, then they don't
    contribute to the distance;
  * if a block from A is similar to a block from B, then mismatched
    parts of these blocks contribute to the distance;
  * if a block from A consists of two blocks in B, then the smaller
    of them contributes to the distance;
  * if a block from A consists of multiple blocks in B, then all of
    them except the biggest one contribute to the distance;
