/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Joiner.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

Joiner::Joiner(int max_dist, float ratio_to_fragment, float gap_ratio):
    max_dist_(max_dist),
    ratio_to_fragment_(ratio_to_fragment),
    gap_ratio_(gap_ratio)
{ }

void Joiner::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("join-max-dist", po::value<int>()->default_value(max_dist()),
     "Max allowed distance when joining fragments")
    ("join-to-fragment", po::value<float>()->default_value(ratio_to_fragment()),
     "Max allowed gap length to fragment length ratio when joining fragments")
    ("join-to-gap", po::value<float>()->default_value(gap_ratio()),
     "Max allowed ratio of gaps' lengths (inside a block) when joining")
   ;
}

void Joiner::apply_options_impl(const po::variables_map& vm) {
    set_max_dist(vm["join-max-dist"].as<int>());
    set_gap_ratio(vm["join-to-gap"].as<float>());
    set_ratio_to_fragment(vm["join-to-fragment"].as<float>());
}

static struct BlockGreater {
    bool operator()(const Block* b1, const Block* b2) const {
        return b1->size() > b2->size();
    }
} block_greater;

static Block* neighbor_block(Block* b, int ori) {
    Block* result = 0;
    Fragment* f = b->front();
    if (f) {
        Fragment* neighbor_f = ori == 1 ? f->next() : f->prev();
        if (neighbor_f) {
            result = neighbor_f->block();
        }
    }
    return result;
}

bool Joiner::can_join(Fragment* one, Fragment* another) {
    return one->seq() == another->seq() && one->ori() == another->ori() &&
           one->is_neighbor(*another);
}

int Joiner::can_join(Block* one, Block* another) {
    bool all[3] = {true, false, true};
    for (int ori = 1; ori >= -1; ori -= 2) {
        BOOST_FOREACH (Fragment* f, *one) {
            Fragment* f1 = f->logical_neighbor(ori);
            if (!f1 || f1->block() != another || !Joiner::can_join(f, f1)) {
                all[ori + 1] = false;
                break;
            }
        }
        if (all[ori + 1]) {
            break;
        }
    }
    int result = all[1 + 1] ? 1 : all[-1 + 1] ? -1 : 0;
    BOOST_ASSERT(!(result && !Block::match(one, another)));
    return result;
}

Block* Joiner::join(Block* one, Block* another, int logical_ori) {
    BOOST_ASSERT(Joiner::can_join(one, another) == logical_ori);
    Block* result = new Block();
    std::set<Fragment*> to_delete;
    BOOST_FOREACH (Fragment* f, *one) {
        Fragment* f1 = f->logical_neighbor(logical_ori);
        BOOST_ASSERT(f1);
        result->insert(join(f, f1));
        to_delete.insert(f);
        to_delete.insert(f1);
    }
    BOOST_FOREACH (Fragment* f, to_delete) {
        delete f;
    }
    return result;
}

Fragment* Joiner::join(Fragment* one, Fragment* another) {
    BOOST_ASSERT(Joiner::can_join(one, another));
    if (another->next() == one) {
        std::swap(one, another);
    }
    Fragment* new_fragment = new Fragment(one->seq());
    new_fragment->set_min_pos(std::min(one->min_pos(), another->min_pos()));
    new_fragment->set_max_pos(std::max(one->max_pos(), another->max_pos()));
    new_fragment->set_ori(one->ori());
    if (one->prev()) {
        Fragment::connect(one->prev(), new_fragment);
    }
    if (another->next()) {
        Fragment::connect(new_fragment, another->next());
    }
    return new_fragment;
}

bool Joiner::can_join_fragments(Fragment* f1, Fragment* f2) const {
    BOOST_ASSERT(Joiner::can_join(f1, f2));
    int dist = f1->dist_to(*f2);
    int min_length = std::min(f1->length(), f2->length());
    BOOST_ASSERT(min_length > 0);
    float ratio = float(dist) / float(min_length);
    return (max_dist_ == -1 || dist <= max_dist_) &&
           (ratio_to_fragment_ < 0 || ratio <= ratio_to_fragment_);
}

bool Joiner::can_join_blocks(Block* b1, Block* b2) const {
    BOOST_ASSERT(Joiner::can_join(b1, b2));
    BOOST_ASSERT(!b1->empty() && !b2->empty());
    Fragment* neighbor_1 = b1->front()->logical_neighbor(1);
    int ori = (neighbor_1 && neighbor_1->block() == b2) ? 1 : -1;
    BOOST_ASSERT(b1->front()->logical_neighbor(ori)->block() == b2);
    int min_gap = -1, max_gap = -1;
    BOOST_FOREACH (Fragment* f1, *b1) {
        Fragment* f2 = f1->logical_neighbor(ori);
        if (!can_join_fragments(f1, f2)) {
            return false;
        }
        int dist = f1->dist_to(*f2);
        min_gap = (min_gap == -1 || dist < min_gap) ? dist : min_gap;
        max_gap = (max_gap == -1 || dist > max_gap) ? dist : max_gap;
    }
    if (gap_ratio_ >= 0 && float(max_gap) / float(min_gap) > gap_ratio_) {
        return false;
    }
    return true;
}

Block* Joiner::try_join(Block* one, Block* another) const {
    Block* result = 0;
    int match_ori = Block::match(one, another);
    if (match_ori == -1) {
        another->inverse();
    }
    if (match_ori) {
        int logical_ori = Joiner::can_join(one, another);
        if (logical_ori && can_join_blocks(one, another)) {
            result = join(one, another, logical_ori);
        }
    }
    return result;
}

bool Joiner::run_impl() const {
    bool result = false;
    std::vector<Block*> bs(block_set()->begin(), block_set()->end());
    std::sort(bs.begin(), bs.end(), block_greater);
    BOOST_FOREACH (Block* block, bs) {
        if (block_set()->has(block)) {
            for (int ori = -1; ori <= 1; ori += 2) {
                while (Block* other_block = neighbor_block(block, ori)) {
                    Block* new_block = try_join(block, other_block);
                    if (new_block) {
                        block_set()->erase(block);
                        block_set()->erase(other_block);
                        block_set()->insert(new_block);
                        block = new_block;
                        result = true;
                    } else {
                        break;
                    }
                }
            }
        }
    }
    return result;
}

}

