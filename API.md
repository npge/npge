TODO: introduction

Words of this document formatted `like this`,
represent names of functions or classes
of the program or inline source code.

NPGe-script language involves data classes (models)
and processors (algorithms).

## Models

Sorry, this section is incomplete.

> [Lua in 15 minutes](http://tylerneylon.com/a/learn-lua)

> If you do not see results of expressions in Lua
> terminal, use function `print`.

There are several data classes used to represent state
of nucleotide pangenome (or, more general, blockset):

 - `Sequence` stores string representing genome sequence
    (e.g., chromosome, plasmid or contig).
    It can represent temporary generated sequence
    (e.g., consensus of a block).
 - `Fragment` points to some fragment of Sequence
    (direct or reverse).
 - `Block` is a collection of fragments.
    A block can store non-aligned fragments
    or aligned fragments.
    In the latter case, each fragment points to
    corresponding instance of AlignmentRow.
 - `AlignmentRow` corresponds to a Fragment and stores
    positions of gaps, representing a row from alignment.
 - `BlockSet` is a collection of Blocks. A BlockSet can be
    a nucleotide pangenome if it satisfies the nucleotide
    pangenome criteria (see file README.md).
    BlockSet points collection of sequences used
    while this BlockSet construction.
 - `BSA` (blockset alignment) is alignment of fragments
    (see screenshot of tool `qnpge`,
    right top part of window).
    If a BlockSet is a pre-pangenome, then
    it can have one or more blockset alignments.
    It is represented in the program as
    alignment. "Letters" of this alignment are
    instances of Fragment, "sequences" of this alignment are
    ordered lists of fragments from one sequence.
    Blockset alignment may include gaps.

![Graphical User Interface of NPG-explorer](http://i.imgur.com/f1LNSSL.png)

Few words about blockset alignment.
If blockset is a pre-pangenome,
then blockset alignment can be constructed and stored in
BSA object.
Blockset alignment is an alignment of fragments themselves,
not their sequences.
In case of pre-pangenome, a sequence is covered
by fragments without intersections.
As in sequence alignments, a row of blockset alignment
represents a sequence or its contiguous part
as a list of fragments with inserted "gap fragments".
In visualization each fragment is replaced with identifier
of the block containing this fragment,
see picture above.
Objective function of blockset alignment awards columns,
in which all fragments are from the same block and
penalizes columns with different blocks and gaps.

Objects of these classes are created with `new` static
method:

```lua
> block_set = BlockSet.new()
> sequence = Sequence.new()
```

Instances of `BSA` can't be created directly from Lua.
See below section about processors.

### Deletion of objects

> This behaviour can be changed in the future
> (I hope I can find a way to make these objects
> garbage collectable).

Objects of `Fragment`, `Block` and `AlignmentRow` require
manual deletion. If you create instance of one of these
classes, you are responsible for deleting it manually
or transferring ownership to other object which in turn
is deleted automatically.
`Fragment` owns corresponding `AlignmentRow`.
`Block` (if it is not `weak`) owns all Fragments it has.
`BlockSet` owns all Blocks it has.
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

By *detaching* we mean unlinking some object
from its owner without deletion of unlinked object.
Example: detaching AlignmentRow from Fragment.

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
Optionally it may include the following value:
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

`Fragment` represents a part of sequence.
`Fragment` is defined by `Sequence`, minimum position
in sequence, maximum position in sequence (0-based) and
orientation (a.k.a. `ori`, can be 1 and -1).

Fragments created on Lua side, must be deleted
manually or attached to other objects (Blocks).

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
(may include gaps ('-') if the fragment has `AlignmentRow`.
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

Method `substr` returns string, representing a subfragment
(as part of a sequence, i.e. withour gaps), it takes
arguments: minimum and maximum positions in fragment
(`minimum <= maximum`).
A negative number passed as a fragment position,
is replaced with `length - abs(value)`.

```lua
> sequence = Sequence.new('ATGC')
> fragment = Fragment.new(sequence, 1, 3)
> fragment:substr(1, -1)
"GC"
```

Method `subfragment` returns a subfragment as new instance
of `Fragment` (make sure to own it properly).
This method takes positions in parent Fragment:
`from` and `to`. If `from > to`, then orientation
of new fragment is opposite to orientation of parent Fragment.

```lua
> sequence = Sequence.new('ATGC')
> fragment = Fragment.new(sequence, 1, 3)
> fragment1 = fragment:subfragment(0, 1)
> fragment1:contents()
"TG"
```

Method `clone` returns new copy of the Fragment.
The copy is not owned by anybody, so make sure
to own it properly. Alignment row is copied as well.

```lua
> sequence = Sequence.new('ATGC')
> fragment = Fragment.new(sequence, 1, 3)
> fragment1 = fragment:clone()
> fragment1:contents()
"TGC"
```

Method `id` returns identifier of a fragment. The
identifier is a string `<seq>_<start>_<stop>`, where
`<seq>` is the sequence name,
`<start>` is the position of fragment begin
in the sequence,
`<stop>` is the position of fragment end
in the sequence.
If `stop < start`, then orientation of the fragment
is `-1`.
If length of fragment is 1 and orientation is `-1`,
then `stop` is replaced with `-1` to distinguish
such a fragment from the fragment with orientation `1`.

```lua
> sequence = Sequence.new('ATGC')
> sequence:set_name('ABC')
> fragment = Fragment.new(sequence, 1, 3)
> fragment:id()
"ABC_1_3"
```

Method `hash` returns hash of underlying part of sequence:

```lua
> sequence = Sequence.new('ATGC')
> fragment = Fragment.new(sequence, 1, 3)
> fragment:hash()
57
```

Method `valid` returns if boundaries of the fragment
are valid (inside Sequence).

Method `has` takes sequence position and return if
it is inside the fragment.

Method `raw_at` takes fragment position and
return nucleotide of this position.
The fragment position can be outside the fragment.
Method `at` is like `raw_at`, but it interprets
negative positions relative to fragment end:

```lua
> sequence = Sequence.new('ATGC')
> fragment = Fragment.new(sequence, 1, 3)
> fragment:at(-1)
"C"
> fragment:raw_at(-1)
"A"
```

Method `common_positions` takes another fragment
and returns number of nucleotides which belong to
block fragments. This method is useful to determine
if two fragment overlap.

```lua
> seq = Sequence.new('ATGC')
> f1 = Fragment.new(seq, 0, 2)
> f2 = Fragment.new(seq, 1, 3)
> f1:common_positions(f2)
2
```

Method `common_fragment` takes another fragment
and returns new instance of `Fragment`, which
occupies nucleotides belonging to both fragments.
Orientation of new fragment is equal to orientation
of the fragment to which method `common_fragment`
is applied.

```lua
> seq = Sequence.new('ATGC')
> f1 = Fragment.new(seq, 0, 2)
> f2 = Fragment.new(seq, 1, 3)
> f3 = f1:common_fragment(f2)
> f3:contents()
"TG"
```

Method `dist_to` takes another fragment and returns
number of nucleotides between them:

```lua
> seq = Sequence.new('ATGC')
> f1 = Fragment.new(seq, 0, 0)
> f2 = Fragment.new(seq, 2, 3)
> f1:dist_to(f2)
2
```

Method `is_subfragment_of` takes another fragment and returns
if the first fragment is inside the second one. This means,
the first fragment has no nucleotides which are not in
the second fragment. If the fragments are equal, this
method returns `true`.

```lua
> seq = Sequence.new('ATGC')
> f1 = Fragment.new(seq, 1, 1)
> f2 = Fragment.new(seq, 1, 2)
> f1:is_subfragment_of(f2)
true
> f2:is_subfragment_of(f1)
false
> f1:is_subfragment_of(f1)
true
```

Method `is_internal_subfragment_of` takes another fragment and
returns if this fragment belongs to another one and
does not share boundaries with it.
If this method returns true, then `is_subfragment_of`
must also be true.

```lua
> seq = Sequence.new('ATGC')
> f1 = Fragment.new(seq, 1, 1)
> f2 = Fragment.new(seq, 1, 2)
> f1:is_internal_subfragment_of(f2)
false
> f2:is_internal_subfragment_of(f1)
false
> f1:is_internal_subfragment_of(f1)
false
> f3 = Fragment.new(seq, 0, 3)
> f2:is_internal_subfragment_of(f3)
true
```

#### Comparison of fragments

Fragments can be compared with `==` and `<`.
Fragments are equal if tuples composed of
the following properties
are equal: minimum and maximum positions,
orientation and sequence (as C pointer).
Fragment `a` is compared to be less than fragment `b`,
if corresponding tuple is less.

#### Methods involving alignment information

Fragment has some methods which use alignment information.
*Block position* is the same as *alignment position*.

- `alignment_at` returns nucleotide by its alignment
    position (like `at` by fragment position, see
    example below)
- `alignment_length` returns length of the alignment row
    (and the alignment itself)
- `set_row` and `row` - access associated `AlignmentRow`
- `detach_row` detaches and returns associated `AlignmentRow`
- `block_pos` takes position in fragment
    and returns corresponding position in block
- `fragment_pos` takes position in block
    and returns corresponding position in fragment
    (if this block position is gap,
    returns nearest position if fragment)
- `frag_to_seq` takes position in fragment
    and returns corresponding position in sequence
- `seq_to_frag` takes position in sequence
    and returns corresponding position in fragment

```lua
> seq = Sequence.new('ATGC')
> f = Fragment.new(seq, 1, 3)
> f:set_row(AlignmentRow.new('T-GC'))
> f:length()
3
> f:alignment_length()
4
> f:block_pos(0) -- T
0
> f:block_pos(1) -- G
2
> f:fragment_pos(0) -- T
0
> f:fragment_pos(1) -- gap => T
0
> f:fragment_pos(2) -- G
1
> f:frag_to_seq(0) -- T
1
> f:seq_to_frag(1) -- T
0
```

### AlignmentRow

AlignmentRow represents string of gaps and non-gaps
(one row of an alignment).
Fragment owns corresponding AlignmentRow.

To create new instance of `AlignmentRow`, use class
function `new`. Optionally you can pass sequence
to the row.
AlignmentRows created on Lua side, must be deleted
manually or attached to other objects (Fragments).

```lua
> row = AlignmentRow.new('T-GC')
```

AlignmentRow doesn't remember actual sequence,
it only remembers at which positions gaps are placed.

To grow right end of AlignmentRow, use method `grow`:

```lua
> row:grow('AA--TT')
```

Another way of adding information to AlignmentRow
is method `bind`. It takes position in fragment
and position in alignment row and "binds" them,
so the fragment position is mapped to the row position.
This method can result in inserting gaps before the position,
if necessary.

To clear the row, use method `clear`:

```lua
> row:clear()
```

To get length of the alignment row, use method `length`:

```lua
> row:length()
10
```

To get alignment position by fragment position,
use method `map_to_alignment`.
To get fragment position by alignment position,
use method `map_to_fragment`. If the position
is a gap, then returns `-1`.
Method `nearest_in_fragment` returns nearest
non-gap fragment position by alignment position.

