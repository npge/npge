/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "OverlaplessUnion.hpp"
#include "FragmentCollection.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

namespace npge {

struct BlockLengthLess {
    bool operator()(Block* a, Block* b) const {
        typedef boost::tuple<int, int, const std::string&> Tie;
        return Tie(b->size(), b->alignment_length(), b->name())
               <
               Tie(a->size(), a->alignment_length(), a->name());
    }
};

static bool check_ou(const OverlaplessUnion* p,
                     std::string& message) {
    bool move = p->opt_value("ou-move").as<bool>();
    bool filter = p->opt_value("ou-filter").as<bool>();
    if (move && filter) {
        message = "--ou-move and --ou-filter are incompatible";
        return false;
    }
    return true;
}

OverlaplessUnion::OverlaplessUnion() {
    add_opt("ou-move", "Move blocks from other to target", false);
    add_opt("ou-filter", "Do not copy good blocks to target, "
            "but remove bad blocks from other", false);
    add_opt_check(boost::bind(check_ou, this, _1));
    declare_bs("target", "Destination of blocks addition");
    declare_bs("other", "Source of blocks addition");
}

void OverlaplessUnion::run_impl() const {
    bool move = opt_value("ou-move").as<bool>();
    bool filter = opt_value("ou-filter").as<bool>();
    BlockSet& t = *block_set();
    BlockSet& o = *other();
    SetFc s2f;
    s2f.add_bs(t);
    Blocks blocks(o.begin(), o.end());
    std::sort(blocks.begin(), blocks.end(), BlockLengthLess());
    BOOST_FOREACH (Block* block, blocks) {
        bool overlaps = s2f.block_has_overlap(block);
        if (overlaps && filter) {
            o.erase(block);
        }
        if (!overlaps && !filter) {
            s2f.add_block(block);
            if (move) {
                o.detach(block);
                t.insert(block);
            } else {
                t.insert(block->clone());
            }
        }
    }
}

const char* OverlaplessUnion::name_impl() const {
    return "Add clones of blocks from other, non-overlapping target";
}

}

