/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <queue>
#include <boost/foreach.hpp>

#include "OverlapsResolver.hpp"
#include "Connector.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

bool OverlapsResolver::overlaps() const {
    BOOST_FOREACH (Block* block, *block_set()) {
        BOOST_FOREACH (Fragment* fragment, *block) {
            for (int ori = -1; ori <= 1; ori += 2) {
                Fragment* neighbor = fragment->neighbor(ori);
                if (neighbor && fragment->common_positions(*neighbor)) {
                    return true;
                }
            }
        }
    }
    return false;
}

struct BlockLess {
    BlockLess(BlockSet* block_set):
        block_set_(block_set)
    { }

    bool operator()(const Block* b1, const Block* b2) const {
        return (block_set_->has(b1) && block_set_->has(b2)) ?
               b1->size() < b2->size() : false;
    }

    BlockSet* block_set_;
};

typedef std::priority_queue<Block*, std::vector<Block*>, BlockLess> BQ;

static void treat_fragments(BlockSet* block_set, BQ& bs,
                            Fragment* x, Fragment* y) {
    Block* x_block = x->block();
    Block* y_block = y->block();
    if (x_block == y_block) {
        x_block->erase(x);
        return;
    }
    Fragment common = x->common_fragment(*y);
    BOOST_ASSERT(common.valid());
    if (*x == common && x->length() == y->length()) {
        BOOST_ASSERT(y_block);
        x->block()->merge(y_block);
        BOOST_ASSERT(y_block->empty());
        block_set->erase(y_block);
    } else {
        if (common == *x) {
            treat_fragments(block_set, bs, y, x);
        } else {
            size_t new_length;
            if (common.begin_pos() == x->begin_pos()) {
                new_length = common.length();
            } else {
                new_length = std::min(abs(x->begin_pos() - common.min_pos()),
                                      abs(x->begin_pos() - common.max_pos()));
            }
            Block* new_block = x->block()->split(new_length);
            BOOST_ASSERT(new_block && !new_block->empty());
            bs.push(new_block);
            block_set->insert(new_block);
        }
    }
}

static bool treat_block(BlockSet* block_set, BQ& bs, Block* block) {
    BOOST_FOREACH (Fragment* f, *block) {
        for (int ori = -1; ori <= 1; ori += 2) {
            Fragment* o_f = f->neighbor(ori);
            if (o_f && f->common_positions(*o_f)) {
                treat_fragments(block_set, bs, f, o_f);
                return true;
            }
        }
    }
    return false;
}

bool OverlapsResolver::run_impl() const {
    bool result = false;
    BQ bs(block_set()->begin(), block_set()->end(),
          BlockLess(block_set().get()));
    while (!bs.empty()) {
        Block* block = bs.top();
        bs.pop();
        while (block_set()->has(block) &&
                treat_block(block_set().get(), bs, block)) {
            result = true;
        }
    }
#ifndef NDEBUG
    BOOST_ASSERT(!overlaps());
    Connector connector;
    connector.apply(block_set());
    BOOST_ASSERT(!overlaps());
#endif
    return result;
}

}

