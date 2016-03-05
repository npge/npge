# Greedy multiple alignment

We introduce our own alignment algorithm which is suitable for
very similar sequences.

**Input**: a list of sequences.

**Output**: multiple alignment.

**Parameters**:

 * (integer) `MISMATCH_CHECK`, minimum number of equal columns
   around a single mismatch (default: 1);

 * (integer) `GAP_CHECK`, minimum number of equal columns around
   a single gap (default: 2);

 * (integer) `ANCHOR`, initial length of anchor used to restart
   alignment (default: 7);

 * (integer) `ANCHOR_FRAME`, maximum position to look for anchor
   (default: 100);

 * (boolean) `RIGHT_UNALIGNED`. If it is set, then the end of the
   input is considered unaligned. It is used when extending a block
   by adding outside fragments: in this case, the beginning is
   always aligned (because the block boundary is aligned), but
   the end can be unaligned because of gaps in some sequences.
   (Default: `false`.)

**Procedures**:

 * moveIdentical
 * alignLeft
 * main

## Procedure moveIdentical

It takes a list of sequences and finds the longest identical
block starting from the beginning of every sequence.
It returns this block and the remaining parts of the sequences.

Example:

```
atagcaat      ata   gcaat
ataccacc   -> ata + ccacc
atagcaggg     ata   gcaggg
```

## Procedure alignLeft

It takes a list of sequences and finds the longest block of high
similarity starting from the beginning of every sequence.
It returns this block and the remaining parts of the sequences.

Steps:

 1. create a list of indices to current positions in input sequences.
    Initialize them with zeros. (Position indices are 0-based.)

 2. If all current positions exist and the pointed letters are equal,
    add current positions to the answer and increment all indices
    and go to step 2.

 3. If all current positions exist and are surrounded by at least
    `MISMATCH_CHECK` identical positions, add current positions to
    the answer and increment all indices and go to step 2.

 4. If `GAP_CHECK` positions to the left of current position are
    identical, try to insert a gap of length 1 to some sequences to get
    `GAP_CHECK` identical positions around the current position.
    The details of this step are discussed below.
    On success, add current positions to the answer
    and increment all indices and go to step 2.

 5. Stop the algorithm.
    The remaining parts of input are considered unaligned.

If some positions do not exist (are outside of the input),
consider them to be identical at steps 3 and 4.

How the gap is found. All letters on the current position are considered.
Given letter `x`, insert a gap before all sequences, which have `x` at
current position. Count the number of identical positions after current
position in the alignment. If it is greater of equal to `GAP_CHECK`,
this choice is a candidate to be the answer. If there is multiple
candidates, choose one with the most number of identical columns.

Examples:

```
atagcaat      atagca   at
ataccacc   -> atacca + cc
atagcaggg     atagca   ggg

ataccaat      atacca   at
ataccacc   -> atacca + cc
ata-caggg     ata-ca   ggg

MISMATCH_CHECK = GAP_CHECK = 2.
```

## Procedure findAnchor

It takes a list of sequences and initial anchor length.

It finds equal fragments of length in range `[GAP_CHECK, ANCHOR]`
which is called anchor.
Longer anchors are preferred. If there are multiple variants,
the variant with all positions equal is preffered. If there is
no variant with all positions equal or there are multiple of them,
the variant with lower values in positions' indices is preferred.

Returns if the anchor was found and on success, returns also
three lists of sequences:

 * parts of sequences before the anchor,
 * anchor (n equal strings),
 * parts of sequences after the anchor.

Steps:

 1. set the current anchor length `L` to `ANCHOR`.
 2. Create a map from word to a list of positions.
    This map stores where a word was seen.
 3. Iterate word beginning from 0 to maximum value for which a word
    of length `L` exists in every sequence and which is less or equal
    to `ANCHOR_FRAME`.
   * Make a list of words at current word beginning `I`.
   * If all words are equal, stop the algorithm: the answer is found.
   * Add pairs (word, `I`) to the map. If a word exists in all sequences,
     stop the algorithm: the answer is found.
 4. Decrement the current anchor length. If it is below `GAP_CHECK`,
    stop the algorithm and return failure. Otherwise go to step 2.

Performance of this step can be improved using a suffix tree.

Examples:

```
Anchor can be found on different positions:

    attattccccccccg -> ------- + attattc + cccccccg
    gggggggattattca    ggggggg   attattc   a-------

Prefer same positions:

    attattcattattcg -> attattc + attattc + g
    gggggggattattca    ggggggg   attattc   a

Length of the anchor can be less than `ANCHOR`:

    cgcgttcattattc -> cgcgttc + atta + ttc
    gggggggattacct    ggggggg   atta   cct

Sometimes no anchor is found:

    attattcattattc -> no anchor
    gggggggcgcgcct

ANCHOR = 7. GAP_CHECK = 2.
```

## Main procedure

 1. Apply `moveIdentical` to the input. If not `RIGHT_UNALIGNED`, then
    apply `moveIdentical` to the complemented form of the remaining parts.
    Add identical part of the result to the answer.

 2. Apply alignLeft to the remaining parts. If not `RIGHT_UNALIGNED`, then
    apply alignLeft to the complemented form of the remaining parts.
    Add aligned part of the result to the answer.

 3. Find anchor in the remaining part using findAnchor. If no anchor was
    found, convert the remaining part to an alignment as-is by adding gaps.
    If an anchor was found, add it to the answer. In this case, also add the
    part before the anchor as-is and reapply the main procedure to the part
    after the anchor.

## Improvements

There are some improvements of the algorithm above:

 1. While adding gaps to the remaining part to convert it to alignment,
    try adding them from the beginning and from the end. Choose the
    alignment with higher identity.

 2. Refine alignment. Move every letter in alignment near a gap across
    the gap. Choose the variant with higher identity.

## Implementation

This algorithm was implemented partly in C++ and Lua in library
lua-npge. You can try it in NPGe tool:

```lua
$ npge
> npge = require 'npge'
> npge.alignment.alignRows({'ATAGATA', 'TAGATAC', 'TACATAC'})
{
    "ATAGATA-",
    "-TAGATAC",
    "-TACATAC",
}
```

Warning! Input must contan only letters "A", "T", "G", "C" and "N"
in upper case.
