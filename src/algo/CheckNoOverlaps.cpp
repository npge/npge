/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "CheckNoOverlaps.hpp"
#include "OverlapsResolver.hpp"
#include "Connector.hpp"
#include "Exception.hpp"

#include "BlockSet.hpp"

namespace npge {

class CheckNoOverlapsImpl : public OverlapsResolver {
protected:
    void run_impl() const {
        if (overlaps()) {
            throw Exception("Overlaps detected");
        }
    }
};

CheckNoOverlaps::CheckNoOverlaps() {
    add(new CheckNoOverlapsImpl);
    add(new Connector);
    add(new CheckNoOverlapsImpl);
    declare_bs("target", "Target blockset");
}

const char* CheckNoOverlaps::name_impl() const {
    return "Make sure no blocks overlap other blocks, else throw";
}

}

