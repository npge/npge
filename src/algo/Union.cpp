/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Union.hpp"
#include "Block.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "BlockSet.hpp"

namespace npge {

Union::Union(const BlockSetPtr& source) {
    set_other(source);
    declare_bs("other", "Source of blocks copying");
    declare_bs("target", "Destination of blocks copying");
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
    block_set()->add_sequences(other()->seqs());
}

const char* Union::name_impl() const {
    return "Add clones of blocks from other to this blockset";
}

}

