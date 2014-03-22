/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
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

namespace bloomrepeats {

typedef std::set<Fragment*, FragmentCompare> FragmentsSet;
typedef FragmentCollection<Fragment*, FragmentsSet> S2F;

static void add_b(S2F& s2f, Block* b) {
    s2f.add_block(b);
}

static void add_bs(S2F& s2f, const BlockSet& bs) {
    s2f.add_bs(bs);
}

static bool f_overlaps(const S2F& s2f, Fragment* f) {
    return s2f.has_overlap(f);
}

static bool b_overlaps(const S2F& s2f, Block* b) {
    BOOST_FOREACH (Fragment* f, *b) {
        if (f_overlaps(s2f, f)) {
            return true;
        }
    }
    return false;
}

typedef std::vector<Block*> Blocks;

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
}

void OverlaplessUnion::run_impl() const {
    bool move = opt_value("ou-move").as<bool>();
    bool filter = opt_value("ou-filter").as<bool>();
    BlockSet& t = *block_set();
    BlockSet& o = *other();
    S2F s2f;
    add_bs(s2f, t);
    Blocks blocks(o.begin(), o.end());
    std::sort(blocks.begin(), blocks.end(), BlockLengthLess());
    BOOST_FOREACH (Block* block, blocks) {
        bool overlaps = b_overlaps(s2f, block);
        if (overlaps && filter) {
            o.erase(block);
        }
        if (!overlaps && !filter) {
            add_b(s2f, block);
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

