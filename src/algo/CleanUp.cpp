/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "CleanUp.hpp"
#include "Joiner.hpp"
#include "Connector.hpp"
#include "OverlapsResolver.hpp"
#include "FragmentsExpander.hpp"
#include "BlocksExpander.hpp"
#include "Filter.hpp"

namespace bloomrepeats {

CleanUp::CleanUp() {
    using namespace boost;
    shared_ptr<Filter> filter = make_shared<Filter>();
    add(filter); // filter.set_min_fragment_length(10);
    add(new Connector);
    add(filter);
    shared_ptr<OverlapsResolver> resolver = make_shared<OverlapsResolver>();
    add(resolver);
    add(new Joiner(0)); // FIXME options
    add(filter);
    add(new BlocksExpander);
    add(resolver);
    add(new FragmentsExpander);
    add(new Filter(100)); // FIXME options
    add(new Joiner(1000));
}

}

