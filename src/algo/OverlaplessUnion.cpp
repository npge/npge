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

#include "OverlaplessUnion.hpp"
#include "FragmentCollection.hpp"
#include "Union.hpp"
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
        return a->alignment_length() < b->alignment_length();
    }
};

bool OverlaplessUnion::run_impl() const {
    bool result = false;
    BlockSet& t = *block_set();
    const BlockSet& o = *other();
    S2F s2f;
    add_bs(s2f, t);
    Blocks blocks(o.begin(), o.end());
    std::sort(blocks.begin(), blocks.end(), BlockLengthLess());
    BOOST_FOREACH (Block* block, blocks) {
        if (!b_overlaps(s2f, block)) {
            add_b(s2f, block);
            t.insert(Union::clone_block(block));
            result = true;
        }
    }
    return result;
}

const char* OverlaplessUnion::name_impl() const {
    return "Add clones of blocks from other, non-overlapping target";
}

}

