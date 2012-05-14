/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <map>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

void BlockSet::insert(BlockPtr block) {
#ifndef NDEBUG
    BOOST_FOREACH (BlockPtr b, *this) {
        BOOST_ASSERT(block != b);
    }
#endif
    blocks_.insert(block);
}

void BlockSet::erase(BlockPtr block) {
    blocks_.erase(block);
}

size_t BlockSet::size() const {
    return blocks_.size();
}

bool BlockSet::empty() const {
    return blocks_.empty();
}

bool BlockSet::has(BlockPtr block) const {
    return blocks_.find(block) != blocks_.end();
}

void BlockSet::clear() {
    blocks_.clear();
}

BlockPtr BlockSet::front() const {
    return empty() ? BlockPtr() : *(begin());
}

BlockSet::iterator BlockSet::begin() {
    return blocks_.begin();
}

BlockSet::const_iterator BlockSet::begin() const {
    return blocks_.begin();
}

BlockSet::iterator BlockSet::end() {
    return blocks_.end();
}

BlockSet::const_iterator BlockSet::end() const {
    return blocks_.end();
}

struct FragmentCompare {
    bool operator()(const FragmentPtr& f1, const FragmentPtr& f2) const {
        return f1->min_pos() < f2->min_pos() ||
               (f1->min_pos() == f2->min_pos() &&
                f1->max_pos() < f2->max_pos());
    }
} fragment_compare;

void BlockSet::connect_fragments() {
    typedef std::vector<FragmentPtr> Fs;
    typedef std::map<SequencePtr, Fs> Seq2Fs;
    Seq2Fs seq2fs;
    BOOST_FOREACH (const BlockPtr& block, *this) {
        BOOST_FOREACH (const FragmentPtr& fragment, *block) {
            seq2fs[fragment->seq()].push_back(fragment);
        }
    }
    BOOST_FOREACH (Seq2Fs::value_type& seq_and_fs, seq2fs) {
        Fs& fs = seq_and_fs.second;
        std::sort(fs.begin(), fs.end(), fragment_compare);
        for (int i = 1; i < fs.size(); i++) {
            Fragment::connect(fs[i - 1], fs[i]);
        }
    }
}

std::ostream& operator<<(std::ostream& o, const BlockSet& block_set) {
    BOOST_FOREACH (BlockPtr block, block_set) {
        o << *block << std::endl;
    }
    return o;
}

}

