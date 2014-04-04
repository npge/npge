/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <vector>
#include <boost/foreach.hpp>

#include "SameChr.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

SameChr::SameChr() {
    declare_bs("target", "Target blockset");
}

bool SameChr::same_chr(const Block* block) {
    if (block->empty()) {
        return true;
    }
    std::string chr = block->front()->seq()->chromosome();
    BOOST_ASSERT(!chr.empty()); // warning here is better
    BOOST_FOREACH (Fragment* f, *block) {
        BOOST_ASSERT(f->seq());
        if (f->seq()->chromosome() != chr) {
            return false;
        }
    }
    return true;
}

void SameChr::run_impl() const {
    BlockSet& bs = *block_set();
    const std::vector<Block*> blocks((bs.begin()), bs.end());
    BOOST_FOREACH (Block* b, blocks) {
        if (!same_chr(b)) {
            bs.erase(b);
        }
    }
}

const char* SameChr::name_impl() const {
    return "Filter out blocks fragments located on different chromosomes";
}

}

