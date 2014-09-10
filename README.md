# NPG-explorer, Nucleotide PanGenome explorer

## Instruction

### Download NPG-explorer

[Download](https://bitbucket.org/starius/npg-explorer/downloads)
prebuild static executables for Windows.

> If Wine complains about missing `wldap32.dll`,
> wine-ldap.
> On Debian: `sudo apt-get install libwine-ldap:i386`.

> Sorry, currently no static version for Linux is
> available (and dynamic is useless if your environment
> differs from mine). Ironically, the development and even
> compilation for Windows is driven by Linux workstation.
> If you use Linux, please give Wine a chance or
> compile the program manually (see below).

### Prepare fasta files with genomes

#### Input file: table of genomes

Create file of the form:

```
CP003176 BRUAO chr1 c Brucella abortus A13334 chr 1
CP003177 BRUAO chr2 c Brucella abortus A13334 chr 2
CP003174 BRUCA chr1 c Brucella canis HSK A52141 chr 1
CP003175 BRUCA chr2 c Brucella canis HSK A52141 chr 2
```

Fields are:

- chromosome entry identifier in EMBL or RefSeqN
    (database is automatically detected from identifier);
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
`npge` will create files and sub-folders
in current directory.
You can change location of output files using command
line options. To see all options, add `-h` to a command.
To set path to table file (instead of `genomes.tsv`),
pass option `--table` to commands
GetFasta, GetGenes and Rename.

#### Prepare sequences and genes

Run the following command:

```bash
$ npge Prepare
```

The following files are created by this command:

- `genomes-renamed.fasta` is FASTA file with
    genomes on with a nucleotide pangenome
    is to be built;
- `features.bs` is a block set of genes.
    One gene is represented as one block.

Files `genomes-raw.fasta` and `features.embl` contain
unprocessed data from a database.
You can safely remove them.

### Build nucleotide pangenome

```bash
$ npge Pangenome
```

This command creates file `pangenome.bs`.
The file is in Block Set format.

#### Check nucleotide pangenome (optional)

```bash
$ npge CheckPangenome
```

This command makes sure that nucleotide pangenome
in file `pangenome.bs` satisfies
pangenome criteria.

The command prints if the pangenome is Ok and
may print some comments about the pangenome.

This step is done by post-processing as well,
result is saved to file `check/isgood`.

### Run post-processing of nucleotide pangenome

```bash
$ npge PostProcessing
```

This command produces many files, some of them
are located in sub-folders.

Files `*.bs` contain blocksets,
`*.bi` contain tables of blocks' properties,
`*.bsa` contain block set alignments.

  * `pangenome.bs` pangenome (main output of the program);
  * `pangenome.bsa` blockset alignment. Table file
    representing alignment in which "letters" are fragments
    of pangenome. This file is used by GUI viewer `qnpge`;
  * `pangenome.hash` hash of pangenome;
  * `pangenome.info` summary of pangenome stats
    (average identity of blocks, etc);
  * `split.bi` and `split.bi` blocks, splitted by
    diagnostic positions.
    This file is used by GUI viewer `qnpge`;
  * `low.bi` and `low.bi`
    Subblocks of low identity sliced from pangenome blocks.
    This file is used by GUI viewer `qnpge`;

  * directory `check` files related to pangenome checking;

    * `check/isgood` result of check that the pangenome
      pangenome criteria;
    * `check/consensuses.fasta` consensus sequences of blocks
      passed to BLAST;
    * `hits.blast` output of BLAST all-against-all run
      on consensuses;
    * `all-blast-hits.bs` and `all-blast-hits.bi` all BLAST hits
      as block sets;
    * `good-blast-hits.bs` and `good-blast-hits.bi`
      BLAST hits satisfying pangenome criteria,
      which surpass overlapping blocks from the pangenome.
      If pangenome satisfies the criteria,
      there are no such blocks;
    * `non-internal-hits.bs` and `non-internal-hits.bi`
      BLAST hits satisfying pangenome criteria,
      which don't surpass overlapping blocks from the pangenome;
    * `joined.bs` and `joined.bi` joined subsequent blocks
      satisfying the criteria;

  * `mutations` mutations related files:

    * `mut.tsv` table of all mutations (columns are block,
      fragment, start of mutation, stop of mutation,
      letter(s) in consensus, letter (or gap) in the fragment);
    * `mutseq.fasta` FASTA file with sequences composed of
      columns with mutations (+ 1 columns to right and to left)
      of stem blocks;
    * `mutseq-with-blocks.bs` same as previous + stem blocks
      mapped of these sequences;

  * `trees` tree related files:

    * `distances.tsv` table file with distances between
      fragments of same block;
    * `all_trees.tsv` list of trees of all blocks
      construsted using Neighbour-Joining;
    * `nj-constree.tre` consensus tree constructed from
      Neighbour-Joining trees of individual blocks;
    * `nj.branch` list and weight of branches constructed
      from Neighbour-Joining trees of individual blocks;
    * `upgma-constree.tre` consensus tree constructed from
      UPGMA trees of individual blocks;
    * `upgma.branch` list and weight of branches constructed
      from UPGMA trees of individual blocks;

> **How to view `.tre` files using FigTree**:
> open a file with FigTree, set branch label to
> "Diagnostic positions" in
> pop-up window, go to "Branch Labels" section of left menu,
> enable the section's checkbox. Abstract distances
> between nodes are shown under branches. To show number
> of diagnostic positions between corresponding clades,
> select "Diagnostic positions" in drop-down list.

### View results in graphical user interface

```bash
$ qnpge
```

![Graphical User Interface of NPG-explorer](http://i.imgur.com/jJHZEse.png)

This command uses `pangenome.bs` and some of files
created by PostProcessing.

The program window is splitted to 3 parts:

 - top left is the table of genomes;
 - top right is the blockset alignment (may be absent
   if there is no file `pangenome.bsa`);
 - bottom is alignment of selected block.

Columns of blocks table:

 - number of fragments in a block;
 - length of a block;
 - identity of a block;
 - percentage of 'G' and 'C' nucleotides in the block;
 - number of genes overlapping with a block;
 - number of parts splitted from a block;
 - number of low similarity regions.

You can filter blocks by block name, gene name or their
sequence using input located up to block table.
To hide blocks of one fragment, clock checkbox
"only blocks of >= 2 fragments".
Blocks table can be sorted by any column.

Blockset alignment table shows alignment of fragment
on genomes.
Chromosome can be selected using drop-down list located
up to blockset table.
Each sequence is represented as a row of blockset table.
Name of a sequence and its orientation against the alignment
is written in first column.
Fragments of a sequence are represented by cells of
blockset table.
Fragments of one block are coloured similarly.
Orientation of a fragment against the alignment
is indicated by '<' and '>'.

When you navigate in blocks table and blockset alignment,
the alignment of the corresponding block is shown in bottom
part of the program.
Fragment name is shown left to alignment itself.
Background colors in alignment correspond to nucleotide types.
Name of the selected gene is shown
in read-only input located up to the alignment.
You can disable genes representation completely
by unchecking the checkbox "show genes".
Genes are coloured with foreground color white.
Genes on reverse chain (relatively to the fragment
orientation) are marked with underscore.
Overlapping genes are coloured with purple.
Start codons are coloured with black,
stop codons are coloured with gray.
Consensus of the block is shown up to the alignment.
Identical columns without gaps are coloured with black,
identical columns with gaps are coloured with gray,
non-identical columns are white.
Columns numbers are shown up to consensus.
Columns numbers of low similarity regions are
coloured with red.

You can use arrows keys to navigate through the alignment.
Corresponding fragment is selected in blockset alignment.
Use keys "Home" and "End" to go to first and last columns
of the alignment respectively.
If you "go away" from the alignment, the program
switches to corresponding block.
You go to next gene boundary if you press `Ctrl + Arrow Right`
or `Ctrl + Arrow Left`.
You go to next low similarity region
if you press `Shift + Arrow Right`
or `Shift + Arrow Left`.

To change order of sequences in blockset alignment
and block alignment, select some rows (you can use `Ctrl` to
select multiple rows) and press `Ctrl + Arrow Up`
or `Ctrl + Arrow Down`.

## Requirements of a good pangenome

 - no overlapping blocks;
 - sequences are covered entirely by blocks
    (including 1-fragment blocks);
 - alignment is defined for each block of >= 2 fragments;
 - length of any fragment except minor blocks and 1-fragment
    blocks is greater or equal to `MIN_LENGTH`;
 - identity of any but minor block is greater or equal
    to `MIN_IDENTITY`;
 - identity of first and last `MIN_LENGTH` columns of
    any but minor block is greater or equal to `MIN_IDENTITY`;
 - first and last columns of blocks do not contain gaps
    or dangling letters;
 - blast run on consensuses finds no blocks which satisfy
    above criteria and surpass overlapping blocks
    from pangenome;
 - no subsequent blocks can be joined so that resulting
    block satisfies above criteria.

## Build and Install

Main executables are command line tool src/tool/npge
(or src/tool/npge.exe) and GUI tool src/gui/qnpge
(src/gui/qnpge.exe).

To change compiled-in default settings,
run `ccmake .` in build directory.

To generate config file, run `npge -g npge.conf` and
change generated file `npge.conf`.

### Requirements

 - C++ compiler (C++11 is not needed);
 - Build system: make, cmake;
 - Boost library (tested with version 1.42);
 - Lua 5.1 or 5.2;
 - LuaBind library.

Optional:

 - Qt library (tested with version 4.8) for GUI;
 - Readline and NCurses libraries for advanced Lua terminal;
 - Doxygen to build documentation;
 - Markdown builder (e.g. Pandoc) to make README.html;
 - a viewer for trees in Newick format (FigTree, MEGA),
    to view phylogenetic trees of genomes.

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

Windows executables are cross-compiled from Linux
using MinGW cross-compiler.

See file `windows/build-requirements`.

Install [MXE](http://mxe.cc) requirements.

Run `windows/build.sh`. It builds needed windows libraries,
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

To create ZIP file and Installation Wizard for Windows,
install pandoc, wine and nsis and run `windows/package.sh`.

## Model

Sorry, this section is incomplete.

> [Lua in 15 minutes](http://tylerneylon.com/a/learn-lua)

> If you do not see results of expressions in Lua
> terminal, use function `print`.

There are several classes used to represent state
of nucleotide pangenome (or, mo general, block set):

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
    a nucleotide pangenome is it satisfies the nucleotide
    pangenome criteria.
    `BlockSet` keeps collection of
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

## Changelog

 - **Version 0.1.1**. Bugfix:

    - fix running MergeUnique in Pangenome,
    - MergeUnique's crash on unique fragments
        located between two fragments of one block was fixed,
    - fixed several bugs when search in GUI by sequence,
        containing 'N's,
    - add icon npge.ico, use it in GUI,
    - Processor memorizes global Meta in its consctuctor,
    - fix Lua function arg_value,
    - Stem is now exact in Info,
    - NEWICK representation of trees is now multiline,
    - set default number of threads to number of cores,
    - windows: build.sh uses current source dir,
    - windows: script package.sh and Installation Wizard,
    - typos in README fixed.

 - **Version 0.1.0**. Features developed in Summer 2014
 were incorporated in ver. 0.1.0.
 It was published prior to ECCB'14 event.

## Links

 - [homepage](http://mouse.belozersky.msu.ru/tools/npge.html),
 - [source](https://bitbucket.org/starius/npg-explorer)
 is hosted on BitBucket,
 - [downloads](https://bitbucket.org/starius/npg-explorer/downloads),
 - [report a bug](https://bitbucket.org/starius/npg-explorer/issues/new).

