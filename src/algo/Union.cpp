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
    Fragment* f1 = new Fragment(*f);
    if (f->row()) {
        f1->set_row(f->row()->clone());
    }
    return f1;
}

Block* Union::clone_block(Block* source) {
    Block* result = new Block(source->name());
    BOOST_FOREACH (Fragment* f, *source) {
        Fragment* f1 = clone_fragment(f);
        result->insert(f1);
    }
    return result;
}

BlockSetPtr Union::clone_block_set(BlockSetPtr block_set) {
    BlockSetPtr result = new_bs();
    Union cloner(block_set);
    cloner.apply(result);
    return result;
}

bool Union::run_impl() const {
    BOOST_FOREACH (Block* block, *other()) {
        block_set()->insert(clone_block(block));
    }
    return !other()->empty();
}

}

