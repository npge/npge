/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "BlockFinder.hpp"
#include "FragmentFinder.hpp"
#include "OverlapFinder.hpp"
#include "BlockSet.hpp"

namespace npge {

BlockFinder::BlockFinder() {
    add_opt("pattern", "Sequence searched for",
            std::string(""));
    ff_ = new FragmentFinder;
    ff_->set_parent(this);
    of_ = new OverlapFinder;
    of_->set_parent(this);
    declare_bs("other", "where to look for");
    declare_bs("target", "where to look for");
}

void BlockFinder::run_impl() const {
    BlockSetPtr pattern_blocks = new_bs();
    pattern_blocks->add_sequences(other()->seqs());
    ff_->set_bs("target", pattern_blocks);
    ff_->run();
    of_->set_bs("bank", other());
    of_->set_bs("pattern", pattern_blocks);
    of_->set_bs("hits", block_set());
    of_->run();
    // clear
    ff_->remove_bs("target");
    of_->remove_bs("bank");
    of_->remove_bs("pattern");
    of_->remove_bs("hits");
}

const char* BlockFinder::name_impl() const {
    return "Locate block by sequence of its fragments";
}

}

