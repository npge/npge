/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "CheckNoOverlaps.hpp"
#include "FragmentCollection.hpp"
#include "Exception.hpp"

#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

namespace npge {

class CheckNoOverlapsImpl : public Processor {
protected:
    void run_impl() const {
        if (overlaps()) {
            throw Exception("Overlaps detected");
        }
    }

    bool overlaps() const {
        TimeIncrementer ti(this);
        SetFc fc;
        BOOST_FOREACH (Block* block, *block_set()) {
            if (fc.block_has_overlap(block)) {
                return true;
            }
            fc.add_block(block);
        }
        return false;
    }
};

CheckNoOverlaps::CheckNoOverlaps() {
    add(new CheckNoOverlapsImpl);
    declare_bs("target", "Target blockset");
}

const char* CheckNoOverlaps::name_impl() const {
    return "Make sure no blocks overlap other blocks, "
           "else throw";
}

}

