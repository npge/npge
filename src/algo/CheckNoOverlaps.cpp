/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "CheckNoOverlaps.hpp"
#include "OverlapsResolver.hpp"
#include "Connector.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

class CheckNoOverlapsImpl : public OverlapsResolver {
protected:
    bool run_impl() const {
        BOOST_ASSERT(!overlaps());
        return false;
    }
};

CheckNoOverlaps::CheckNoOverlaps() {
#ifndef NDEBUG
    add(new CheckNoOverlapsImpl);
    add(new Connector);
    add(new CheckNoOverlapsImpl);
#endif
}

}

