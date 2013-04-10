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
#include "Union.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

struct FragmentPtrLess {
    bool operator()(Fragment* a, Fragment* b) const {
        return *a < *b;
    }
};

typedef std::set<Fragment*, FragmentPtrLess> FragmentsSet;
typedef std::map<Sequence*, FragmentsSet> S2F;

static void add_f(S2F& s2f, Fragment* f) {
    Sequence* seq = f->seq();
    s2f[seq].insert(f);
}

static void add_b(S2F& s2f, Block* b) {
    BOOST_FOREACH (Fragment* f, *b) {
        add_f(s2f, f);
    }
}

static void add_bs(S2F& s2f, const BlockSet& bs) {
    BOOST_FOREACH (Block* b, bs) {
        add_b(s2f, b);
    }
}

static bool f_overlaps(const S2F& s2f, Fragment* f) {
    Sequence* seq = f->seq();
    S2F::const_iterator it = s2f.find(seq);
    if (it == s2f.end()) {
        return false;
    }
    const FragmentsSet& fragments = it->second;
    BOOST_ASSERT(!fragments.empty());
    FragmentsSet::const_iterator i2 = fragments.lower_bound(f);
    if (i2 != fragments.end() && (*i2)->common_positions(*f)) {
        return true;
    } else if (i2 != fragments.begin()) {
        i2--;
        if ((*i2)->common_positions(*f)) {
            return true;
        }
    }
    return false;
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

}

