/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Union.hpp"
#include "Block.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

Union::Union(const BlockSetPtr& source) {
    set_other(source);
}

Fragment* Union::clone_fragment(Fragment* f) {
    return f->clone();
}

Block* Union::clone_block(Block* source) {
    return source->clone();
}

BlockSetPtr Union::clone_block_set(BlockSetPtr block_set) {
    return block_set->clone();
}

void Union::run_impl() const {
    other()->copy(*block_set());
}

}

