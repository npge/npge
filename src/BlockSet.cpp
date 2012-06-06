/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <queue>
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

static struct FragmentCompare {
    bool operator()(const FragmentPtr& f1, const FragmentPtr& f2) const {
        return *f1 < *f2;
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

void BlockSet::filter(int min_fragment_length, int min_block_size) {
    std::vector<BlockPtr> block_set_copy(begin(), end());
    BOOST_FOREACH (const BlockPtr& block, block_set_copy) {
        block->filter(min_fragment_length);
        if (block->size() < min_block_size) {
            erase(block);
        }
    }
}

static struct BlockGreater {
    bool operator()(const BlockPtr& b1, const BlockPtr& b2) const {
        return b1->size() > b2->size();
    }
} block_greater;

static BlockPtr neighbour_block(const BlockPtr& b, int ori) {
    BlockPtr result;
    FragmentPtr f = b->front();
    if (f) {
        FragmentPtr neighbour_f = ori == 1 ? f->next() : f->prev();
        if (neighbour_f) {
            result = neighbour_f->block();
        }
    }
    return result;
}

void BlockSet::join() {
    std::vector<BlockPtr> bs(begin(), end());
    std::sort(bs.begin(), bs.end(), block_greater);
    BOOST_FOREACH (BlockPtr block, bs) {
        if (has(block)) {
            for (int ori = -1; ori <= 1; ori += 2) {
                while (BlockPtr other_block = neighbour_block(block, ori)) {
                    BlockPtr new_block = Block::try_join(block, other_block);
                    if (new_block) {
                        erase(block);
                        erase(other_block);
                        insert(new_block);
                        block = new_block;
                    }
                }
            }
        }
    }
}

void BlockSet::expand_blocks(PairAligner* aligner, int batch,
                             int ori, bool overlap) {
    std::vector<BlockPtr> bs(begin(), end());
    std::sort(bs.begin(), bs.end(), block_greater);
    BOOST_FOREACH (BlockPtr block, bs) {
        block->expand(aligner, batch, ori, overlap);
    }
}

bool BlockSet::intersections() const {
    BOOST_FOREACH (BlockPtr block, *this) {
        BOOST_FOREACH (FragmentPtr fragment, *block) {
            for (int ori = -1; ori <= 1; ori += 2) {
                FragmentPtr neighbour = fragment->neighbour(ori);
                if (neighbour && fragment->common_positions(*neighbour)) {
                    return true;
                }
            }
        }
    }
    return false;
}

typedef std::map<Fragment, FragmentPtr> F2F;

static void get_middle(const FragmentPtr& fr, const FragmentPtr& intersection,
                       F2F& f2f) {
    FragmentDiff diff = fr->diff_to(*intersection);
    BOOST_FOREACH (FragmentPtr f, *fr->block()) {
        FragmentPtr new_f = boost::make_shared<Fragment>();
        new_f->apply_coords(*f);
        new_f->patch(diff);
        if (!new_f->valid()) {
            continue;
        }
        new_f->set_ori(1); // same for all fragments of new block
        if (f2f.find(*new_f) == f2f.end()) {
            new_f->find_place(f);
            f2f[*new_f] = new_f;
        }
    }
}

static bool inside(FragmentPtr small, BlockPtr large) {
    for (int ori = -1; ori <= 1; ori += 2) {
        FragmentPtr neighbour = small;
        while (neighbour && neighbour->common_positions(*small)) {
            if (neighbour->block() == large) {
                return true;
            }
            neighbour = neighbour->neighbour(ori);
        }
    }
    return false;
}

static bool inside(BlockPtr small, BlockPtr large) {
    BOOST_FOREACH (FragmentPtr fragment, *small) {
        if (!inside(fragment, large)) {
            return false;
        }
    }
    return true;
}

void BlockSet::patch_block(const BlockPtr& block, const FragmentDiff& diff) {
    block->patch(diff);
    block->find_place();
    block->filter(/* min_length */ 1);
    if (block->empty()) {
        erase(block);
    }
#ifndef NDEBUG
    BOOST_FOREACH (FragmentPtr fragment, *block) {
        BOOST_ASSERT(fragment->valid());
    }
#endif
}

static BlockPtr split_block(const FragmentPtr& f, const FragmentPtr& common) {
    FragmentDiff diff = f->diff_to(*common);
    F2F f2f;
    BOOST_FOREACH (FragmentPtr fragment, *f->block()) {
        Fragment middle;
        middle.apply_coords(*fragment);
        middle.patch(diff);
        if (!middle.valid() || !middle.is_internal_subfragment_of(*fragment)) {
            continue;
        }
        Fragment left_f;
        left_f.apply_coords(middle);
        left_f.set_min_pos(fragment->min_pos());
        FragmentPtr right_f = boost::make_shared<Fragment>();
        right_f->apply_coords(middle);
        right_f->set_max_pos(fragment->max_pos());
        fragment->apply_coords(left_f);
        fragment->find_place();
        right_f->find_place(fragment);
        f2f[*right_f] = right_f;
    }
    BlockPtr right = Block::create_new();
    BOOST_FOREACH (F2F::value_type& f_and_ptr, f2f) {
        FragmentPtr& ptr = f_and_ptr.second;
        right->insert(ptr);
    }
    return right;
}

BlockPtr BlockSet::treat_two(const FragmentPtr& x, const FragmentPtr& y,
                             int min_intersection) {
    FragmentPtr intersection = x->common_fragment(*y);
    BOOST_ASSERT(intersection);
    if (x->is_internal_subfragment_of(*y)) {
        return split_block(y, intersection);
    } else if (y->is_internal_subfragment_of(*x)) {
        return split_block(x, intersection);
    }
    BlockPtr result;
    FragmentPtr small_f = x->block()->size() < y->block()->size() ? x : y;
    BlockPtr small = small_f->block();
    FragmentPtr large_f = small_f == x ? y : x;
    BlockPtr large = large_f->block();
    if (intersection->length() >= min_intersection && !inside(small, large)) {
        F2F f2f;
        get_middle(x, intersection, f2f);
        get_middle(y, intersection, f2f);
        result = Block::create_new();
        BOOST_FOREACH (F2F::value_type& f_and_ptr, f2f) {
            FragmentPtr& ptr = f_and_ptr.second;
            result->insert(ptr);
        }
        patch_block(large, large_f->exclusion_diff(*intersection));
    }
    patch_block(small, small_f->exclusion_diff(*intersection));
    BOOST_ASSERT(!x || !y || !x->common_positions(*y));
    return result;
}

static struct BlockLess {
    bool operator()(const BlockPtr& b1, const BlockPtr& b2) const {
        return b1->size() < b2->size();
    }
} block_less;

void BlockSet::resolve_intersections(int min_intersection) {
    std::priority_queue<BlockPtr, std::vector<BlockPtr>, BlockLess> bs(begin(),
            end(), block_less);
    while (!bs.empty()) {
        BlockPtr block = bs.top();
        bs.pop();
new_block:
        if (has(block)) {
            BOOST_FOREACH (FragmentPtr f, *block) {
                for (int ori = -1; ori <= 1; ori += 2) {
                    FragmentPtr o_f = f->neighbour(ori);
                    if (o_f && f->common_positions(*o_f)) {
                        BlockPtr b = treat_two(f, o_f, min_intersection);
                        if (b) {
                            insert(b);
                            bs.push(b);
                        }
                        goto new_block;
                    }
                }
            }
        }
    }
#ifndef NDEBUG
    BOOST_ASSERT(!intersections());
    connect_fragments();
    BOOST_ASSERT(!intersections());
#endif
}

bool BlockSet::expand_blocks_by_fragments(PairAligner* aligner) {
    bool result = false;
    BOOST_FOREACH (BlockPtr block, *this) {
        BOOST_ASSERT(block);
        result |= block->expand_by_fragments(aligner);
    }
    return result;
}

std::ostream& operator<<(std::ostream& o, const BlockSet& block_set) {
    BOOST_FOREACH (BlockPtr block, block_set) {
        o << *block << std::endl;
    }
    return o;
}

}

