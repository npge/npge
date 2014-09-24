/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "CheckNoOverlaps.hpp"
#include "Connector.hpp"
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
        BOOST_FOREACH (Block* block, *block_set()) {
            BOOST_FOREACH (Fragment* fragment, *block) {
                for (int ori = -1; ori <= 1; ori += 2) {
                    Fragment* neighbor =
                        fragment->neighbor(ori);
                    if (neighbor && fragment->common_positions(
                                *neighbor)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
};

CheckNoOverlaps::CheckNoOverlaps() {
    add(new Connector);
    add(new CheckNoOverlapsImpl);
    declare_bs("target", "Target blockset");
}

const char* CheckNoOverlaps::name_impl() const {
    return "Make sure no blocks overlap other blocks, else throw";
}

}

