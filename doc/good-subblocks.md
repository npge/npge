# NPG-explorer

## How good subblocks are detected in a bad block

### Parameters

 * `min_length`
 * `min_identity`

These two options are used while constructing the pangenome.
Meaning of `min_length` is a minimum length of a pangenome
block. Typical value of `min_length` is 100. Meaning of
`min_identity` is a percentage of identical columns in a
pangenome block. Typical value of `min_identity` is 0.9.

### Input

A block, which doesn't satisfy the criteria (bad block).

### Outputs

A list of subblocks which satisfy the criteria (good
subblocks).

### Algorithm

1. Each column of a block is categorized to ge "good" or "bad".
   A column is considered "good", if it contains no gaps and
   is represented by the same letter in all fragments.

2. If the block starts or ends with bad columns, exclude
   them from subsequent steps.

3. Iterate all positions of a slice of length `min_length`
   in the block. Calculate identity for a subblock, defined
   by a slice position. If identity of a subblock is greater
   than or equal to `min_identity` then the corresponding
   slice is considered "good".

4. Group continuous good frames into larges frames.

5. Create set of selected good subblocks (empty).

6. Find good joined slice of maximum length.
   If its length is greater than or equal to `min_length`,
   convert it into a subblock and add this subblock to
   the output set of good subblocks.
   Otherwise stop. The algorithm has finished

7. Exclude positions of the slice selected from all
   remaining good frames.

8. Go to step 6.

Step 3 can be done in `O(N)`, where N is the block's length.
Create counter which keeps current number of good columns.
Iterate first `min_length` columns, add number of good columns
to the counter. Check it's value: if it is greater than or
equal to `min_identity * min_length`, then add slice
`[1; min_length]` to set of good frames. "Move" the slice
forward: add value of next column to the counter and subtract
value of the column leaving the slice. Apply the same check
to the slice. Repeat this step `N - min_length` times.

Selection of good frame of maximum length (step 6) can be
optimised using binary heap, binary tree or similar data
structure.

### Example

Parameters:
 - `min_length` = 3
 - `min_identity` = 0.5

Input:
```
TAATGC
TATT-C
```

The algorithm's steps are as follows:

```
TAATGC
TATT-C

++-+-+   Detect good and bad columns:

+++
 +++     Detect good and bad frames of length 3:
  ---
   +++

++++     Join continuous good frames, create a subblock:
```
