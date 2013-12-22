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
#include <boost/algorithm/string/join.hpp>

#include "OneByOne.hpp"
#include "DeConSeq.hpp"
#include "Filter.hpp"
#include "Align.hpp"
#include "Union.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "hit.hpp"

namespace bloomrepeats {

typedef std::vector<Block*> Blocks;

struct BlockSizeLess {
    bool operator()(Block* a, Block* b) const {
        return b->size() < a->size();
    }
};

struct OneByOne::Impl {
    Filter* filter_;
    Align* align_;

    Impl(OneByOne* obo) {
        filter_ = new Filter;
        filter_->set_parent(obo);
        align_ = new Align;
        align_->set_parent(obo);
    }
};

OneByOne::OneByOne() {
    impl_ = new Impl(this);
}

OneByOne::~OneByOne() {
    delete impl_;
}

static bool has_overlap(S2F& s2f, Block* block) {
    BOOST_FOREACH(Fragment* f, *block) {
        if (s2f.has_overlap(f)) {
            return true;
        }
    }
    return false;
}

static void insert_or_delete(Block* block, BlockSet& target, S2F& s2f) {
    if (!block->empty()) {
        BOOST_ASSERT(!has_overlap(s2f, block));
        target.insert(block);
        s2f.add_block(block);
    } else {
        delete block;
    }
}

struct LeftRight {
    Block* left;
    Block* right;

    LeftRight() {
        left = new Block;
        right = new Block;
    }
};

typedef std::map<Block*, LeftRight> Overlap2LR;

static void split_fragment(Overlap2LR& o2lr, Fragment* fragment,
                           Fragment* overlap_fragment) {
    Block* overlap_block = overlap_fragment->block();
    BOOST_ASSERT(overlap_block);
    LeftRight& lr = o2lr[overlap_block];
    std::vector<size_t> positions;
    positions.push_back(overlap_fragment->min_pos());
    positions.push_back(overlap_fragment->max_pos());
    positions.push_back(fragment->min_pos());
    positions.push_back(fragment->max_pos());
    std::sort(positions.begin(), positions.end());
    if (positions[0] != positions[1]) {
        Fragment* min_fragment = new Fragment(overlap_fragment->seq(),
                                              positions[0], positions[1] - 1,
                                              overlap_fragment->ori());
        Block* container = (overlap_fragment->ori() == 1) ?
                           lr.left : lr.right;
        container->insert(min_fragment);
    }
    if (positions[2] != positions[3]) {
        Fragment* max_fragment = new Fragment(overlap_fragment->seq(),
                                              positions[2] + 1, positions[3],
                                              overlap_fragment->ori());
        Block* container = (overlap_fragment->ori() == -1) ?
                           lr.left : lr.right;
        container->insert(max_fragment);
    }
}

static void split_block(S2F& s2f, Block* hit, BlockSet& target) {
    Overlap2LR o2lr;
    Block* hit_clone = Union::clone_block(hit);
    typedef std::map<Fragment*, Fragment*> F2F;
    F2F overlap_to_hit;
    std::set<Block*> orig_blocks;
    std::set<std::string> names;
    BOOST_FOREACH (Fragment* fragment, *hit) {
        std::vector<Fragment*> overlap_fragments;
        s2f.find_overlap_fragments(overlap_fragments, fragment);
        BOOST_ASSERT(overlap_fragments.size() == 1);
        Fragment* overlap_fragment = overlap_fragments[0];
        orig_blocks.insert(overlap_fragment->block());
        names.insert(overlap_fragment->block()->name());
        F2F::iterator hit_fr = overlap_to_hit.find(overlap_fragment);
        if (hit_fr != overlap_to_hit.end()) {
            hit_fr->second = 0; // >=2 hit fragments in one old fragment.
            // maybe, in-block repeat
        } else {
            overlap_to_hit[overlap_fragment] = fragment;
        }
    }
    BOOST_FOREACH (F2F::value_type& overlap_and_hit, overlap_to_hit) {
        Fragment* overlap_fragment = overlap_and_hit.first;
        Fragment* hit = overlap_and_hit.second;
        if (hit) {
            split_fragment(o2lr, hit, overlap_fragment);
        }
    }
    BOOST_FOREACH (Block* block, orig_blocks) {
        s2f.remove_block(block);
    }
    BOOST_FOREACH (Overlap2LR::value_type& overlap_and_lr, o2lr) {
        Block* orig_block = overlap_and_lr.first;
        Block* left = overlap_and_lr.second.left;
        Block* right = overlap_and_lr.second.right;
        left->set_name(orig_block->name() + "_left");
        right->set_name(orig_block->name() + "_right");
        insert_or_delete(left, target, s2f);
        insert_or_delete(right, target, s2f);
    }
    BOOST_FOREACH (Block* block, orig_blocks) {
        target.erase(block);
    }
    hit_clone->set_name(boost::algorithm::join(names, "_"));
    insert_or_delete(hit_clone, target, s2f);
}

bool OneByOne::run_impl() const {
    bool result = false;
    BlockSet& t = *block_set();
    BlockSet& o = *other();
    Blocks hits_blocks(o.begin(), o.end());
    std::sort(hits_blocks.begin(), hits_blocks.end(), BlockSizeLess());
    S2F s2f;
    s2f.add_bs(t);
    Align* align = impl_->align_;
    Filter* filter = impl_->filter_;
    BOOST_FOREACH (Block* hit, hits_blocks) {
        filter->filter_block(hit);
        if (!filter->is_good_block(hit)) {
            continue;
        }
        if (!is_internal_hit(s2f, hit)) {
            continue;
        }
        align->apply_to_block(hit);
        if (!filter->is_good_block(hit)) {
            continue;
        }
        if (!is_internal_hit(s2f, hit)) {
            continue;
        }
        if (has_self_overlaps(hit)) {
            fix_self_overlaps(hit);
            if (!filter->is_good_block(hit)) {
                continue;
            }
            align->apply_to_block(hit);
            if (!filter->is_good_block(hit)) {
                continue;
            }
            if (!is_internal_hit(s2f, hit)) {
                continue;
            }
        }
        split_block(s2f, hit, t);
        result = true;
    }
    return result;
}

const char* OneByOne::name_impl() const {
    return "Split blockset according to hits, processed one by one";
}

}
