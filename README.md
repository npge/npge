[![Build Status](https://travis-ci.org/npge/npge.png?branch=master)](https://travis-ci.org/npge/npge)
[![License](https://img.shields.io/badge/License-GPL2,3-brightgreen.png)](LICENSE)

## Instruction

### Download NPG-explorer

[Download](https://github.com/npge/npge/releases)
prebuild static executables for Windows and Linux.

BLAST and other dependencies (except Qt 4
in Linux version) are included.

In the following instructions, replace `x.y.z` with
the version of NPG-explorer you use.

For Windows, download and run file `npge_x.y.z_win32.exe`
as administrator.

For Linux 32 bit, download and unpack file
`npge_x.y.z_lin32.tar.gz`
(using command `tar -xf npge_x.y.z_lin32.tar.gz`).

For Linux 64 bit, download and unpack file
`npge_x.y.z_lin64.tar.gz`
(using command `tar -xf npge_x.y.z_lin64.tar.gz`).

Windows version adds itself to PATH, so you can use
commands `npge` and `qnpge` in command line.
In Linux, you need to add the unpacked directory
`npge-x.y.z` to PATH. If you use `bash`,
open `~/.bashrc` in your favorite text editor
and add the following line to the end:
`export PATH=$PATH:/path/to/npge-x.y.z`.
Do not forget to replace `/path/to/npge-x.y.z`
with the actual path.

### Prepare a directory for the results

All files with data and output files have fixed names.
This is why it is recomended to make separate
working directory for each task.

### Prepare files with genomes

#### Input file: table of genomes

Create file of the form:

```
all:embl:CP003176 BRUAO chr1 c Brucella abortus A13334 chr 1
all:embl:CP003177 BRUAO chr2 c Brucella abortus A13334 chr 2
fasta:refseqn:NC_016778.1 BRUCA chr1 c Brucella canis HSK A52141 chr 1
features:embl:CP003174 BRUCA chr1 c Brucella canis HSK A52141 chr 1
fasta:file:BRUCA.fasta BRUCA chr1 c Brucella canis HSK A52141 chr 1
features:file:BRUCA.fasta BRUCA chr1 c Brucella canis HSK A52141 chr 1
fasta:file:base.fasta[CP002459] BRUMM chr1 c Brucella melitensis
features:file:base.embl[CP002459] BRUMM chr1 c Brucella melitensis
```

Fields are:

- chromosome entry identifier; it is composed from
  record type ('fasta', 'features' (annotation) or 'all'),
  source (database name or 'file')
  and identifier in that database or file path;
  it is possible to get annotations and fasta files
  from different sources;
  use `file.fasta[name_of_sequence]` as an identifier
  to copy specific sequence from file;
- short name for the genome chosen by user,
    this name is used in output data;
    this name must not include spaces or underscores
    (more precisely, it must not look like a fragment
    name: `AAA_NNN_NNN`);
- chromosome name (e.g., 'chr1', 'chr2'),
- chromosome circularity ('c' for circular and 'l'
    for linear),
    and arbitrary description (not used by the program).

Annotations in the following formats are parsed
by the program: GenBank, EMBL.

String `CP003175 BRUCA chr2 c` corresponds to EMBL entry
`CP003175` which is represented by short genome name
`BRUCA`, chromosome name `chr2` and is circular.

You can use contigs instead of chromosomes,
if genome is not fully assembled.
Set circularity to 'l' in this case.

Such a table for 17 genomes of Brucella can be found in
file examples/Brucella/17genomes.tsv.
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
- `genes/features.bs` is a blockset of genes.
    One gene is represented as one block.

Files `genomes-raw.fasta` and `features.embl` contain
unprocessed data from a database.
They are not used by following steps.
You can safely remove them.

#### Examine prepared sequences

Run the following command:

```bash
$ npge Examine
```

The following files are created by this command in
directory `examine`:

- `genomes-info.tsv` table of genome lengths;
    make sure genome lengths are about the same,
    otherwise do not expect NPG-explorer to build
    good pangenome;
- `draft.bs` draft pangenome;
- `identity_recommended.txt` text file with recommended value
    of parameter `MIN_IDENTITY` (see below how to change
    this parameter in configuration file).

This step is needed to gather some information
about input genomes. This information can be used
on next step (configuration).

### Set values of global options

To change values of global options, make file `npge.conf`
using command `npge -g npge.conf`, then edit this file
to change values of global options.
File `npge.conf` contains default values compiled into
the program. Sometimes they have to be changed to
improve results.

Configuration file looks like this:

```
MIN_IDENTITY = Decimal('0.9')
MIN_LENGTH = 100
```

Decimal values are specified using the syntax above.
Accuracy of decimal values is 4 digits after the point.

The program applies following configuration files (if exist):

 - `npge.conf` in current directory;
 - `npge.conf` in the program's directory (directory where
 the executable lives);
 - `/etc/npge.conf`
 - `~/.npge.conf`
 - reads environmental variables named like known options;
 - `npge.conf` in current directory (again);
 - applies command line options if passed.

### Build nucleotide pangenome

```bash
$ npge MakePangenome
```

This command creates file `pangenome/pangenome.bs`.
The file is in BlockSet format.

#### Check nucleotide pangenome (optional)

```bash
$ npge CheckPangenome
```

This command makes sure that nucleotide pangenome
in file `pangenome/pangenome.bs` satisfies
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
`*.ba` contain blockset alignments.
Files `*.bi` contain numbers of occurrences of a block
in a genome. Each genome adds one column to the table.
To add similar columns with numbers of occurrences of
a block in a sequence, add option `--info-count-seqs=1`
to command `npge PostProcessing`.

  * `pangenome/pangenome.bs` pangenome
    (main output of the program);
  * `pangenome/pangenome.ba` blockset alignment. Table file
    representing alignment in which "letters" are fragments
    of pangenome. This file is used by GUI viewer `qnpge`;
  * `global-blocks/blocks.bs` global blocks (blockset).
    Global blocks are joined collinear S-blocks;
  * `global-blocks/blocks.ba` global blocks
    (blockset alignment).
    This file is used by GUI viewer `qnpge`;
  * `global-blocks/blocks.gbi` properties of global
    blocks as a table;
  * `pangenome/pangenome.hash` hash of pangenome;
  * `pangenome/pangenome.info` summary of pangenome stats
    (average identity of blocks, etc);
  * `extra-blocks/split.bs` and `extra-blocks/split.bi`
    blocks, splitted by diagnostic positions.
    This file is used by GUI viewer `qnpge`;
  * `extra-blocks/low.bs` and `extra-blocks/low.bi`
    Subblocks of low identity sliced from pangenome blocks.
    This file is used by GUI viewer `qnpge`;

  * directory `check` files related to pangenome checking;

    * `isgood` result of check that the pangenome
      pangenome criteria;
    * `hits.blast` output of BLAST all-against-all run
      on consensuses;
    * `all-blast-hits.bs` and `all-blast-hits.bi` all BLAST hits
      as blocksets;
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
      of stable blocks;
    * `mutseq-with-blocks.bs` same as previous + stable blocks
      mapped of these sequences;
    * `consensuses.fasta` consensus sequences of blocks;

  * `genes` analysis of genes:

    * `genes/good.bs` groups of genes matching each other
      exactly accourding to S-blocks;
    * `genes/good-upstreams.bs` upstream sequences of
      good genes (`genes/good.bs`);

  * `trees` tree related files:

    * `nj-global-tree.tre` global tree constructed using
      Neighbour-Joining applied to genome distances;

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

![Graphical User Interface of NPG-explorer](http://i.imgur.com/f1LNSSL.png)

This command uses `pangenome/pangenome.bs` and some of files
created by PostProcessing.

The program window is splitted to 3 parts:

 - top left is the table of genomes;
 - top right is the blockset alignment (may be absent
   if there is no file `pangenome/pangenome.ba`);
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
By default, pattern matching is wildcard.
`^` before the pattern and `$` after the pattern
correspond to name start/end (as in regular expressions).
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
    or dangling letters (few letters followed by long gaps);
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
 - Lua (version 5.1 or 5.2) or LuaJIT;
 - LuaBind library;
 - BLAST legacy or BLAST plus (run-time requirement).

Optional:

 - Qt library (tested with version 4.8) for GUI;
 - Readline and NCurses libraries for advanced Lua terminal;
 - Doxygen to build documentation;
 - Markdown builder (e.g. Pandoc) to make README.html;
 - a viewer for trees in Newick format (FigTree, MEGA),
    to view phylogenetic trees of genomes.

### Linux

> To build static Linux package in fresh Debian Wheezy,
> install curl and sudo and run
> `curl -L http://tiny.cc/npge-build-linux | sh`

Install build requirements (on Debian):

```bash
% ./linux/requirements.sh
```

Build the program as static executables
(Qt 4 is not static!):

```bash
$ ./linux/build.sh
```

The program is built in the directory `npge-build-linux`.

Create distribution `.tar.gz` file:
go into `npge-build-linux` and run:

```bash
$ ./linux/package.sh
```

How to build manually:

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Pass argument `-DNPGE_STATIC_LINUX:BOOL=1` to after `cmake`
to get static executables (on Debian, Qt 4 is not static!).

Build README.html:

```bash
$ pandoc -s README.html README.md
```

Run tests:

```bash
$ make test
```

### Windows

> To build static Windows packages in fresh Debian Wheezy,
> install curl and sudo and run
> `curl -L http://tiny.cc/npge-build-windows | sh`

Windows executables are cross-compiled from Linux
using MinGW cross-compiler.

Install build requirements (on Debian):

```bash
% ./windows/requirements.sh
```

Build the program as static executables:

```bash
$ ./windows/build.sh
```

The program is built in the directories
`npge-build-windows32` and `npge-build-windows64`
(first contains 32-bit version and
second contains 64-bit version).

Create ZIP file and Installation Wizard for Windows,
go into build directory (`npge-build-windows32` or
 `npge-build-windows64`) and run:

```bash
$ ./windows/package.sh
```

## Changelog

 - **Version 0.3.1**. Bug fix.
    - fix crashing in processor MergeUnique.

 - **Version 0.3.0**. Improvements.
    - non-compatibles changes of formats,
    - new syntax of of script language added (see API.md),
    - Fragment now has not `prev` and `next`
        (use `FragmentCollection` instead),
    - some processors were renamed,
    - documentation was updated and improved,
    - Windows 64bit build,
    - LuaJIT by default,
    - blockset alignment is build in global blocks,
    - GUI: search by sequence.

 - **Version 0.2.5**. Cleanup, bugfix:
    - fix MergeUnique fail if sequence has only 2 fragment,
    - fix MakePangenome's imperfections on joining,
    - prevent BLAST failures on long N's
    (split before blasting),
    - move FragmentCollection to model/,
    - add new features to FragmentCollection,
    - fix bug in FragmentCollection.remove,
    - start moving code from fragment.prev/next to
    FragmentCollection (Rest, Joiner, bsa, MergeUnique),
    - remove FragmentDiff, fragment.split, fragment.exclude,
    fragment.shift, OneByOne,
    - `Sequence.circular()` doesn't throw (linear by default),
    - `Sequence.set_name()` throws if name includes space
    or underscore,
    - fix memory leak in a test.

 - **Version 0.2.1**. Bugfix:
    - Windows Uninstaller removes all shortcuts from Start,
    - mycoplasma example,
    - fix app path discovery if UPX'ed,
    - wrap cmd with quotes if it contains spaces,
    - BlastFinder: do not crash if input is empty,
    - GetData can copy specific records from local files.

 - **Version 0.2.0**. Cleanup:
    - TrySmth: fix error (subblocks),
    - Pangenome: make sure result blocks can not be joined,
    - remove deprecated and unused code and options,
    - change format of config files,
    - read config `npge.conf` located in application dir,
    - group options to sections and hide some of them,
    - change format of `genomes.tsv` (no automagic,
    can read from file),
    - import other features, not only genes, update GUI,
    - GUI/search: not only wildcard search,
    - GUI: ignore missing additional input files,
    - print "Running...", ".. done",
    - less technical output of Pangenome,
    - rename :cin => :stdin, :cout => :stdout, :cerr => :stderr,
    - remove filter options `spreading` and `max_gaps`.

 - **Version 0.1.4**. Bugfix:

    - clarifications in README,
    - fix parsing errors in `read_block_set` and `AddGenes`,
    - fix possible test failure,
    - fix possible terminal freeze,
    - support alternative Lua implementation LuaJIT,
    - Windows Wizard creates more shortcuts,
    - throw exception on file openning error.

 - **Version 0.1.3**. Bugfix:

    - configure to use blast plus for static build,
    - add blast to requirements,
    - move from hg and bitbucket to git, github,
    - update npge.ico,
    - one-line build commands.

 - **Version 0.1.2**. Bugfix:

    - fix "Pangenome" (some blocks were not joined),
    - fix freeze in `thread_pool`,
    - do not use CURL,
    - fix many build issues,
    - all-steps.sh for Windows and Linux,
    - run tests before packaging,
    - fix all tests to run them from wine,
    - Linux static build + package.

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
 - [source](https://github.com/npge/npge)
 is hosted on GitHub,
 - [downloads](https://github.com/npge/npge/releases)
 - [report a bug](https://github.com/npge/npge/issues/new)
 (GitHub sign-in is needed).

This work was presented at ECCB'14 conference:
[abstract](http://www.eccb14.org/95) and
[poster](http://f1000.com/posters/browse/summary/1096831).

Corresponding author: Boris Nagaev, email: bnagaev@gmail.com

