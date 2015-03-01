/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "SelfOverlapsResolver.hpp"
#include "hit.hpp"

namespace npge {

SelfOverlapsResolver::SelfOverlapsResolver() {
    declare_bs("target", "Target blockset");
}

void SelfOverlapsResolver::process_block_impl(Block* block,
        ThreadData*) const {
    fix_self_overlaps(block);
}

const char* SelfOverlapsResolver::name_impl() const {
    return "Resolve self-overlaps";
}

}

