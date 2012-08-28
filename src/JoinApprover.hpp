/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_JOIN_APPROVER_HPP_
#define BR_JOIN_APPROVER_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Utility object, making decision whether blocks/fragments can be merged.
Blocks/fragments must be joinable (Block::can_join and Fragment::can_join).
*/
class JoinApprover {
public:
    /** Constructor.
    \param max_dist Max allowed \ref Fragment::dist_to "distance".
        Value -1 means that this limitation is not applied.
    */
    JoinApprover(int max_dist = -1);

    /** Return if two fragments can be joined.
    The fragments must be joinable (Fragment::can_join).
    Blocks are not taken into account here.
    */
    bool can_join_fragments(FragmentPtr f1, FragmentPtr f2);

    /** Return if two blocks can be joined.
    The blocks must be joinable (Block::can_join) and non empty.
    Both fragments and blocks are taken into account here.
    */
    bool can_join_blocks(BlockPtr b1, BlockPtr b2);

private:
    int max_dist_;
};

}

#endif

