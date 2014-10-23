/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_BSA_ALGO_HPP_
#define NPGE_BSA_ALGO_HPP_

#include <iosfwd>

#include "block_set_alignment.hpp"
#include "global.hpp"

namespace npge {

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
               const BSA& first, const BSA& second, int genomes);

/** Produce alignment from sub-alignments.
Align sub-alignments one by one.
*/
void bsa_make_aln(BSA& aln, const BSAs& parts, int genomes);

/** Produce alignment from map of trivial rows */
void bsa_make_aln(BSA& aln, const BSA& rows, int genomes);

/** Produce alignment from map of trivial rows using tree
Tree leaf nodes should return sequence names.
*/
void bsa_make_aln_by_tree(BSA& aln, const BSA& rows,
                          const TreeNode* tree, int genomes);

/** Remove pure gap columns from alignment */
void bsa_remove_pure_gaps(BSA& aln);

/** Improve alignment by moving fragments left and right */
void bsa_move_fragments(BSA& aln);

/** Split mixed columns with gaps to several columns */
void bsa_unwind(BSA& aln);

/** Move columns if possible to get neighbours nearby */
void bsa_move_columns(BSA& aln);

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

/** Choose orientation of sequences by majority, move start.
If all sequences are circular, start is moved so that sum
of indices of boundaries is minimum. Fragments must be connected.
*/
void bsa_orient(BSA& bsa);

/** Replace all fragments from non-stem blocks (exact) with gaps */
void bsa_filter_exact_stem(BSA& bsa, int genomes);

/** Replace all fragments from short blocks with gaps */
void bsa_filter_long(BSA& bsa, int min_length);

/** Print blockset alignment.
\param out Output stream.
\param aln Blockset alignment.
\param name Name of blockset alignment.
    It is printed in the beginning of each line.
*/
void bsa_print(std::ostream& out, const BSA& aln,
               const std::string& name,
               bool orientation = true);

/** Print one row with conservative blocks */
void bsa_print_conservative(std::ostream& out, const BSA& aln,
                            const std::string& name);

/** Input blockset alignment */
void bsa_input(BlockSet& bs, std::istream& in);

}

#endif

