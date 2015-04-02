# NPG-explorer

## How low identity subblocks are created

### Parameters

 * `min_length`
 * `min_identity`

These two options are used while constructing the pangenome.
Meaning of `min_length` is a minimum length of a pangenome
block. Typical value of `min_length` is 100. Meaning of
`min_identity` is a percentage of identical columns in a
pangenome block. Typical value of `min_identity` is 0.9.

### Input

A block of a pangenome.

### Outputs

A list of low identity subblocks. A low identity subblock is a
block, sliced from the original block. It corresponds to a
range of columns in the alignment of the original block.

### Algorithm

Each column of a block is categorized to ge "good" or "bad".
A column is considered "good", if it contains no gaps and
is represented by the same letter in all fragments.

Then continuous groups of columns of same value of the flag
are formed. All groups get a weight. Initial weight of "good"
group is equal to its length. Initial weight of "bad" group
is calculated as follows:
```
weight = length / (1 - min_identity).
```

Next step is reducing small groups. This step is based on
[the algorithm][kchp] proposed by Shellia Guberman,
see [1].

 1. If number of groups is 1, then the algorithm stops.
 2. A group with minimum weight is selected.
    If there are more than one group with minimum weight,
    then arbitrary group is selected.
 3. If weight of the group selected is more or equal to
    `min_length` the algorithm stops.
 4. Otherwise this group is joined with its neighbour groups.
    They are replaced with new group with flag value equal to
    the value of flag of neighbour groups and weight value
    equal to sum of weights of joined groups.
 5. Go to step 1.

When the algorithm stops, all "bad" groups are selected.
A "bad" group is converted to low identity subblock.

### References

 1. Губерман Ш.А.
    Неформальный анализ данных в геологии и геофизике. -
    М.: Недра, 1987.

[kchp]: http://www.integro.ru/system/ots/guberman/guberman.htm
