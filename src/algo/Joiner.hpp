/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_JOIN_APPROVER_HPP_
#define NPGE_JOIN_APPROVER_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "Decimal.hpp"

namespace npge {

class MetaAligner;

/** Join subsequent blocks.
Blocks/fragments must be joinable (Block::can_join and Fragment::can_join).

\ref Block::weak() "Weak" blocks can't be joined.
*/
class Joiner : public Processor {
public:
    /** Constructor */
    Joiner();

    /** Return if these fragments can be joined (simple check).
    Fragments can be joined if they share the same sequence and ori
    and \ref is_neighbor "are neighbors".
    */
    bool can_join(Fragment* one, Fragment* another) const;

    /** Return if these blocks can be joined (simple check).
     - 1, if fragments of 'one' should preceed fragments from 'another';
     - -1, if fragments of 'another' should preceed fragments from 'one';
     - 0, if blocks can't be joined.

    Blocks can be joined, if they match and all the fragments from
    the first block has an unique neighbor with the same ori
    from the second block.
    */
    int can_join(Block* one, Block* another) const;

    /** Merge fragments and return new larger fragment */
    Fragment* join(Fragment* one, Fragment* another) const;

    /** Return joined blocks, if these two blocks can be joined */
    Block* join_blocks(Block* one, Block* another,
                       int logical_ori) const;

    /** Return if two fragments can be joined.
    The fragments must be joinable (Fragment::can_join) and non empty.
    Blocks are not taken into account here.
    */
    bool can_join_fragments(Fragment* f1, Fragment* f2) const;

    /** Return if two blocks can be joined.
    The blocks must be joinable (Block::can_join) and non empty.
    Both fragments and blocks are taken into account here.
    */
    bool can_join_blocks(Block* b1, Block* b2) const;

    /** Try to join, return empty pointer if failed */
    Block* try_join(Block* one, Block* another) const;

protected:
    /** Apply the action */
    void run_impl() const;

    const char* name_impl() const;

private:
    void build_alignment(Strings& rows,
                         const Fragments& fragments,
                         const Block* another,
                         int logical_ori) const;
    Block* neighbor_block(Block* b, int ori) const;
    MetaAligner* aligner_;
};

}

#endif

