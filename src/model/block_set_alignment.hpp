/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_SET_ALIGNMENT_HPP_
#define BR_BLOCK_SET_ALIGNMENT_HPP_

#include <iosfwd>
#include <map>
#include <vector>

#include "global.hpp"

namespace bloomrepeats {

/** Array of fragments, 0 for gap */
typedef std::vector<Fragment*> Fragments;

/** One row of block set alignment */
struct BSRow {
    /** Constructor */
    BSRow();

    /** Ori of row in relation to sequence */
    int ori;

    /** Array of fragments, 0 for gap */
    Fragments fragments;
    // TODO circular genomes
    // TODO start of row in arbitrary part of sequence
};

/** Block set alignment */
typedef std::map<Sequence*, BSRow> BSA;

/** Vector of block set alignment */
typedef std::vector<BSA> BSAs;

/** Return length of block set alignment */
int bsa_length(const BSA& bsa);

/** Create blocks set alignment row of the sequence.
Output BSA is not an alignment.
If rows is empty, all sequences of the chromosome will be added.
If rows is not empty, only sequences existing in rows will be added.
*/
void bsa_make_rows(BSA& rows, const BlockSet& bs);

/** Inverse alignment */
void bsa_inverse(BSA& aln);

/** Create blocks set alignment row of the sequence */
void bsa_align(BSA& both, int& score,
               const BSA& first, const BSA& second);

/** Produce alignment from sub-alignments.
Align sub-alignments one by one.
*/
void bsa_make_aln(BSA& aln, const BSAs& parts);

/** Produce alignment from map of trivial rows */
void bsa_make_aln(BSA& aln, const BSA& rows);

/** Produce alignment from map of trivial rows using tree
Tree leaf nodes should return sequence names.
*/
void bsa_make_aln_by_tree(BSA& aln, const BSA& rows,
                          const TreeNode* tree);

/** Produce tree
Ownership is transferred to caller.
Tree leaf nodes return sequence names.
Tree is built using N-J.
Distance between two sequences is calculated as number of
common blocks divided to number of blocks.

\warning Do not delete rows, until you use resulting tree.
*/
TreeNode* bsa_make_tree(const BSA& rows);

/** Convert tree of genomes to tree of sequences.
Ownership is transferred to caller.
Input tree has genome or sequence names as leaf names.
Input tree has sequence names as leaf names.
If two sequences share same genome name, exception is thrown.

\warning Do not delete rows, until you use resulting tree.
*/
TreeNode* bsa_convert_tree(const BSA& rows, const TreeNode* tree);

/** Print block set alignment.
\param out Output stream.
\param aln Block set alignment.
\param name Name of block set alignment.
    It is printed in the beginning of each line.
\param blocks Print block names, else fragments.
*/
void bsa_print(std::ostream& out, const BSA& aln,
               const std::string& name,
               bool blocks = true);

/** Input block set alignment */
void bsa_input(BlockSet& bs, std::istream& in);

}

#endif

