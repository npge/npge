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

static struct BlockCompare {
    bool operator()(const BlockPtr& b1, const BlockPtr& b2) const {
        return b1->size() > b2->size();
    }
} block_compare;

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

void BlockSet::merge() {
    std::vector<BlockPtr> bs(begin(), end());
    std::sort(bs.begin(), bs.end(), block_compare);
    BOOST_FOREACH (BlockPtr block, bs) {
        if (has(block)) {
            for (int ori = -1; ori <= 1; ori += 2) {
                while (BlockPtr other_block = neighbour_block(block, ori)) {
                    BlockPtr new_block = Block::try_merge(block, other_block);
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
    std::sort(bs.begin(), bs.end(), block_compare);
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
        new_f->set_ori(1); // same for all fragments of new block
        if (f2f.find(*new_f) == f2f.end()) {
            if (f->next()) {
                Fragment::connect(new_f, f->next());
            }
            Fragment::connect(f, new_f);
            new_f->find_place();
            f->find_place();
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
}

BlockPtr BlockSet::treat_two(const FragmentPtr& x, const FragmentPtr& y,
                             int min_intersection) {
    FragmentPtr intersection = x->common_fragment(*y);
    BOOST_ASSERT(intersection);
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
    return result;
}

void BlockSet::resolve_intersections(int min_intersection) {
    std::vector<BlockPtr> bs(begin(), end());
    std::sort(bs.begin(), bs.end(), block_compare);
    BOOST_FOREACH (BlockPtr block, bs) {
        if (has(block)) {
new_block:
            BOOST_FOREACH (FragmentPtr f, *block) {
                for (int ori = -1; ori <= 1; ori += 2) {
                    FragmentPtr o_f = f->neighbour(ori);
                    if (o_f && f->common_positions(*o_f)) {
                        BlockPtr block = treat_two(f, o_f, min_intersection);
                        if (block) {
                            insert(block);
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

std::ostream& operator<<(std::ostream& o, const BlockSet& block_set) {
    BOOST_FOREACH (BlockPtr block, block_set) {
        o << *block << std::endl;
    }
    return o;
}

}

