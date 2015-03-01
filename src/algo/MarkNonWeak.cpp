/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "MarkNonWeak.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"

namespace npge {

MarkNonWeak::MarkNonWeak() {
    declare_bs("target", "Target blockset");
}

void MarkNonWeak::run_impl() const {
    BOOST_FOREACH (Block* block, *block_set()) {
        if (block->weak()) {
            block->set_weak(false);
        }
    }
}

const char* MarkNonWeak::name_impl() const {
    return "Mark all blocks non-weak";
}

}

