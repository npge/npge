# NPG-explorer

## Instruction

### Prerare fasta files with genomes

1. Create file of the form:

    ```
    CP003176 BRUAO chr1 c Brucella abortus A13334 chr 1
    CP003177 BRUAO chr2 c Brucella abortus A13334 chr 2
    CP003174 BRUCA chr1 c Brucella canis HSK A52141 chr 1
    CP003175 BRUCA chr2 c Brucella canis HSK A52141 chr 2
    ```

    Fields are:

    - chromosome entry identifier in EMBL or RefSeqN
        (database is automaticaly detected from identifier);
    - short name for the genome chosen by user,
        this name is used in output data;
    - chromosome name (e.g., 'chr1', 'chr2'),
    - chromosome circularity ('c' for circular and 'l'
        for linear),
        and arbitrary description (not used by the program).

    String `CP003175 BRUCA chr2 c` corresponds to EMBL entry
    `CP003175` which is represented by short genome name
    `BRUCA`, chromosome name `chr2` and is circular.

    You can use contigs instead of chromosomes,
    if genome is not fully assembled.
    Set circularity to 'l' in this case.

    Such a table for 17 genomes of Brucella can be found in
    file brucella/17genomes.tsv.
    This example is used below.

    Create empty directory and create file `genomes.tsv`
    with the table of genomes to be used to build pangenome.
    `npge` will create files and subfolders
    in current directory.
    You can change location of output files using command
    line options. To see all options, add `-h` to a command.
    To set path to table file (instead of `genomes.tsv`),
    pass option `--table` to commands
    GetFasta, GetGenes and Rename.

1. Prepare sequences and genes:

    ```bash
    $ npge Prepare
    ```

    The following files are created by this command:

    - `genomes-renamed.fasta` is FASTA file with
        genomes on with a nucleotide pangenome
        is to be built;
    - `features.bs` is a [block set](#blockset) of genes.
        One gene is represented as one block.

    Files `genomes-raw.fasta` and `features.embl` contain
    unprocessed data from a database.
    You can safely remove them.

### Build nucleotide pangenome

```bash
$ npge Pangenome
```

This command creates file `pangenome.bs`.
The file is in [Block Set](#blockset) format.

#### Check nucleotide pangenome (optional)

```bash
$ npge CheckPangenome
```

This command makes sure that nucleotide pangenome
in file `pangenome.bs` satisfies
[pangenome criteria](#criteria).

The command prints if the pangenome is Ok and
may print some comments about the pangenome.

### Run post-processing of nucleotide pangenome

```bash
$ npge PostProcessing
```

This command produces many files, some of them
are located in subfolders.

### View results in graphical user interface

```bash
$ qnpge
```

This command uses `pangenome.bs` and some of files
created by PostProcessing.

## Build and Install

Main executables are command line tool src/tool/npge
(or src/tool/npge.exe) and GUI tool src/gui/qnpge
(src/gui/qnpge.exe).

### Requirements

 - C++ compiler (C++11 is not needed);
 - Build system: make, cmake;
 - Boost library (tested with version 1.42);
 - Lua 5.1 or 5.2;
 - LuaBind library;
 - libCURL.

Optional:

 - Qt library (tested with version 4.8) for GUI;
 - Readline and NCurses libraries for advanced Lua terminal;
 - Doxygen to build documentation;
 - Markdown builder (e.g. Pandoc) to make README.html.

### Linux

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Build README.html:

```bash
$ pandoc -s README.html README.md
```

Run tests:

```bash
$ make test
```

### Windows

See file windows-build-requirements.

Install [MXE](http://mxe.cc) requirements.

Run windows-build.sh. It builds needed windows libraries,
downloads npge-explorer (last version) and builds it.
Executables are linked statically.

Download [BLAST+ binaries](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+-ia32-win32.tar.gz).
Extract files makeblastdb.exe and blastn.exe.
Download
[vcomp100.dll](http://drive5.com/usearch/manual/vcomp100.html),
which BLAST+ executables depend on.

To decrease size of executables install
Ultimate Packer for eXecutables (UPX) and run:

```bash
$ strip *.exe
$ upx -9 *.exe
```

## Model

> [Lua in 15 minutes](http://tylerneylon.com/a/learn-lua)

> If you do not see results of expressions in Lua
> terminal, use function `print`.

There are several classes used to represent state
of nucleotide pangenome (or, mo general, block set):

 - `Sequence` stores string representing genomic sequence.
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
    a nucleotide pangenome is it satisfies the nucleotide
    pangeome criterion. `BlockSet` keeps collection of
    sequences used in its `Fragment`s.
 - `BSA` (block set alignment) stores alignment of fragments.
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

<a name="delete"></a>

Objects of `Fragment`, `Block` and `AlignmentRow` require
manual deletion. If you create instance of one of these
classes, you are responsible for deleting it manually
or transfering ownership to other object which in turn
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

`Sequence` stores string representing genomic sequence.
Moreover, it stores name and description of the sequence.
NPGe recognises the following sequence name format:
`GENOME&CHROMOSOME&CIRCULAR` (3 values joint by '&').
`GENOME` is short name of a genome.
`CHROMOSOME` is short name of a chromosome.
`CIRCULAR` can be 'c' (circular) or 'l' (linear).
This information is used by some algorithms of NPGe.
Example of sequence name: `BRUO2&chr1&c` (chromosome 1 of
Brucella ovis ATCC 25840, which is circular).

Description of a sequence can be any string.
Optionally it can include the following value:
`ac=XXXX`. NPGe recognises this entry as Accession Number
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
> seq:substr(2, 2, -1)
"CA"
```

Note that first argument of `seq:substr(2, 2, -1)`. Because the
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
> seq:hash(2, 2, -1)
"3"
> make_hash("CA")
"3"
```

To filter out any characted but 'A', 'T', 'G', 'C' and 'N'
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

`Fragment`s created on Lua side, must be [deleted
manually](#delete).

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

