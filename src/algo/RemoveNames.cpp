/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "RemoveNames.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Sequence.hpp"

namespace npge {

RemoveNames::RemoveNames() {
    add_opt("remove-blocks-names", "Remove blocks names", true);
    add_opt("remove-seqs-names", "Remove seqences names", true);
    declare_bs("target", "Target blockset");
}

void RemoveNames::run_impl() const {
    if (opt_value("remove-blocks-names").as<bool>()) {
        BOOST_FOREACH (Block* block, *block_set()) {
            block->set_name("");
        }
    }
    if (opt_value("remove-seqs-names").as<bool>()) {
        BOOST_FOREACH (SequencePtr seq, block_set()->seqs()) {
            seq->set_name("");
        }
    }
}

const char* RemoveNames::name_impl() const {
    return "Remove all blocks and/or sequences names";
}

}

