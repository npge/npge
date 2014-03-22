/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "CheckNoOverlaps.hpp"
#include "OverlapsResolver.hpp"
#include "Connector.hpp"
#include "Exception.hpp"

#include "BlockSet.hpp"

namespace bloomrepeats {

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
}

}

