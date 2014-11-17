## Model

Sorry, this section is incomplete.

> [Lua in 15 minutes](http://tylerneylon.com/a/learn-lua)

> If you do not see results of expressions in Lua
> terminal, use function `print`.

There are several classes used to represent state
of nucleotide pangenome (or, more general, blockset):

 - `Sequence` stores string representing genome sequence.
 - `Fragment` points to some fragment of `Sequence`
    (direct or reverse).
 - `Block` is a collection of `Fragment`s. `Fragment`s
    inside block
    can be aligned with each other, in this case alignment
    information is stored in instances of `AlignmentRow`.
 - `AlignmentRow` corresponds to a `Fragment` and stores string
    consisting of gaps and non-gaps, representing a row
    from alignment.
 - `BlockSet` is a collection of `Block`s. A `BlockSet` can be
    a nucleotide pangenome if it satisfies the nucleotide
    pangenome criteria.
    `BlockSet` keeps collection of
    sequences used in its `Fragment`s.
 - `BSA` (blockset alignment) stores alignment of fragments.
    One or more `BSA`s can be stored in `BlockSet`,
    eachof them is accessible by string key.

Objects of classes listed above are long-living.
This means in particular that they can persist across several
`Processor`s (see bellow).
Objects of long-living classes are created with `new` static
method:

```lua
> block_set = BlockSet.new()
> sequence = Sequence.new()
```

Instances of `BSA` can't be created directly, so the only
way to get instance of `BSA` in Lua script is through
methods of `BlockSet`.

### Deletion of objects

Objects of `Fragment`, `Block` and `AlignmentRow` require
manual deletion. If you create instance of one of these
classes, you are responsible for deleting it manually
or transferring ownership to other object which in turn
is deleted automatically.
`Fragment` owns corresponding `AlignmentRow`.
`Block` (if it is not `weak`) owns all `Fragment`s it has.
`BlockSet` owns all `Block`s it has.
To delete an object manually, use `AlignmentRow.delete`,
`Fragment.delete` or `Block.delete`.

```lua
> fragment = Fragment.new()
> Fragment.delete(fragment)
-- do not use fragment after this line
-- unless you want program crash
```

If you create temporary object, you can create a deleter,
which owns this object. When deleter object is garbage
collected, it deletes corresponding object.

```lua
> f2 = fragment:clone()
> df = Fragment.deleter(f2)
-- df deletes f2
```

If you do not delete objects of `AlignmentRow`, `Fragment`
or `Block`, this results in memory leak.
If you use an object after deletion, this results in
program crash.

### Sequence

`Sequence` stores string representing genome sequence.
Moreover, it stores name and description of the sequence.
NPGe recognizes the following sequence name format:
`GENOME&CHROMOSOME&CIRCULAR` (3 values joint by '&').
`GENOME` is short name of a genome.
`CHROMOSOME` is short name of a chromosome.
`CIRCULAR` can be 'c' (circular) or 'l' (linear).
This information is used by some algorithms of NPGe.
Example of sequence name: `BRUO2&chr1&c` (chromosome 1 of
Brucella ovis ATCC 25840, which is circular).

Description of a sequence can be any string.
Optionally it can include the following value:
`ac=XXXX`. NPGe recognizes this entry as Accession Number
of the sequence in databank.

Create new sequence and set its name and description:

```lua
> seq = Sequence.new()
> seq:set_name("TEST&chr1&c")
> seq:set_description("ac=ABC123")
```

Get name, description, genome, chromosome, circularity
and accession number of the sequence:

```lua
> seq:name()
"TEST&chr1&c"
> seq:description()
"ac=ABC123"
> seq:genome()
"TEST"
> seq:chromosome()
"chr1"
> seq:circular()
true
> seq:ac()
"ABC123"
```

`Sequence` stores string consisting of letters
'A', 'T', 'G', 'C' and 'N'. 'N' is used for unknown letter.
To append new part to a Sequence, use method `push_back`:

```lua
> seq:push_back("AT")
> seq:push_back("GC")
```

Create `Sequence` from string:

```lua
> sequence = Sequence.new('ATGC')
```

To get contents of a `Sequence`, use method `contents`:

```lua
> seq:contents()
"ATGC"
```

To get size of a `Sequence`, use method `size`:

```lua
> seq:size()
4
```

If you convert `Sequence` to string (using function
`tostring` or just by entering its name), you get
FASTA representation of this sequence:

```lua
> tostring(seq)
[[
> TEST&chr1&c ac=ABC123
ATGC
]]
```

To get a letter from a sequence by its index (0-based),
use method `char_at`:

```lua
> seq:char_at(0)
"A"
```

To get a substring from a sequence, use method `substr`.
Arguments are start position, length and orientation
(1 is direct orientation, -1 is reversed (complement)
orientation).

```lua
> seq:substr(1, 2, 1)
"TG"
> seq:substr(2, 1, -1)
"CA"
```

Note that first argument of `seq:substr(2, 1, -1)`. Because the
substring is reversed, start position is 2, not 1.

To hash substring of sequence, use method `hash`.
Arguments are the same of `substr`.
Free function `make_hash` can be applied to a string
to get its hash.

```lua
> seq:hash(1, 2, 1)
"9"
> make_hash("TG")
"9"
> seq:hash(2, 1, -1)
"3"
> make_hash("CA")
"3"
```

To filter out any character but 'A', 'T', 'G', 'C' and 'N'
from a string, use `Sequence.to_atgcn`:

```lua
> Sequence.to_atgcn("AT-G-C-A#T$N .C")
"ATGCATNC"
```

There are several ways to store string in `Sequence`:

 - `Sequence.ASIS_SEQUENCE` is plain string (as is).
 One letter takes 1 byte.
 - `Sequence.COMPACT_SEQUENCE` takes advantage of
 alphabet ('A', 'T', 'G', 'C', 'N') to reduce memory
 spent by a sequence.
 One letter takes 3 bits.
 - `Sequence.COMPACT_LOW_N_SEQUENCE` takes advantage
 of fact that 'N' is very rare letter in most sequences.
 One letter takes 2 bits + overhead
 depending on number of 'N's.

By default, one of compact sequences is used.
You can specify sequence type directly:

```lua
> seq = Sequence.new(Sequence.ASIS_SEQUENCE)
> seq = Sequence.new('ATGC', Sequence.ASIS_SEQUENCE)
```

#### DummySequence

`DummySequence` is a sequence all letters in which are equal.

```lua
> seq = DummySequence.new('A', 20)
> seq:contents()
"AAAAAAAAAAAAAAAAAAAA"
```

You can skip letter and/or size of `DummySequence` and
set them by methods `set_letter` and `set_size`. Letter
of `DummySequence` can be retrieved using method `letter`.

```lua
> seq = DummySequence.new()
> seq:set_letter('A')
> seq:set_size(20)
> seq:contents()
"AAAAAAAAAAAAAAAAAAAA"
```

#### FragmentSequence

`FragmentSequence` points to content of a `Fragment`.

```lua
> fragment = ...
> fs = FragmentSequence.new(fragment)
> assert(fs:size() == fragment:size())
> assert(fs:contents() == fragment:str())
```

The fragment can be set latter by method `set_fragment`.

### Fragment

`Fragment` represents part of sequence.
`Fragment` is defined by `Sequence`, minimum position
in sequence, maximum position in sequence (0-based) and
orientation (a.k.a. `ori`, can be 1 and -1).

`Fragment`s created on Lua side, must be deleted
manually.

> Note: `Fragment` and `Sequence` do not own each other,
> though `Fragment` points to `Sequence`.
> If you use the fragment, whose sequence gets deleted,
> the program crashes.

Method `Fragment.new` accepts 4 arguments (`seq`, `min_pos`,
`max_pos`, `ori`), all are optional. They can be accessed
by getters and setters named same as arguments:

```lua
-- create a fragment which covers whole sequence
> fragment = Fragment.new(sequence, 0, sequence:size() - 1)
> assert(fragment:length() == sequence:size())
> fragment:set_min_pos(2)
> fragment:set_max_pos(7)
> assert(fragment:length() == 7 - 2 + 1)
> fragment:ori()
1
> fragment:set_ori(-1)
> fragment:ori()
-1
```

Use method `inverse` to change `ori` to the opposite.

```lua
> fragment:ori()
-1
> fragment:inverse()
> fragment:ori()
1
```

Reversed fragment points reverse chain of the sequence,
thus it contains complement sequence.

Method `begin_pos` returns position in the sequence,
where the fragment starts.
Method `last_pos` returns position in the sequence,
where the fragment ends.
Method `end_pos` returns position in the sequence,
next to the position where the fragment ends and not
located in the fragment. `end_pos` can be located even
not in the sequence.

```lua
> f = Fragment.new(s, 10, 20)
> f:begin_pos()
10
> f:last_pos()
20
> f:end_pos()
21
> f.inverse()
> f:begin_pos()
20
> f:last_pos()
10
> f:end_pos()
9
```

To set `begin_pos` and `last_pos` at once, use method
`set_begin_last`. Orientation is changed to satisfy both
begin and last positions.

```lua
> f = Fragment.new(s)
> f:set_begin_last(10, 20)
> f:ori()
1
> f:set_begin_last(20, 10)
> f:ori()
-1
```

Method `seq` returns sequence of the fragment,
method `block` returns block of the fragment.

```lua
> seq = fragment:seq()
> block = fragment:block()
```

Use method `str` to get part of `Sequence` occupied by
the fragment.
Use method `contents` to get string contents of `Fragment`
(can include gaps ('-') if the fragment has `AlignmentRow`.
Method `length` returns length of part of the sequence
occupied by the fragment.
Method `alignment_length` returns length of contents of
the fragment (including gaps).


```lua
> sequence = Sequence.new('ATGC')
> fragment = Fragment.new(sequence, 1, 2)
> fragment:contents()
"TG"
> fragment:str()
"TG"
> fragment:length()
2
> fragment:alignment_length()
2
> fragment:set_row(AlignmentRow.new('T-G'))
> fragment:contents()
"T-G"
> fragment:str()
"TG"
> fragment:length()
2
> fragment:alignment_length()
3
```

