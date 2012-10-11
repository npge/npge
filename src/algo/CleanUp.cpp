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
    filter->set_min_fragment_length(10);
    filter->set_no_options(true);
    add(filter);
    add(new Connector);
    add(filter);
    shared_ptr<OverlapsResolver> resolver = make_shared<OverlapsResolver>();
    add(resolver);
    shared_ptr<Joiner> joiner = make_shared<Joiner>(0);
    joiner->set_no_options(true);
    add(joiner);
    add(filter);
    add(new BlocksExpander);
    add(resolver);
    add(new FragmentsExpander);
    add(new Filter);
    add(new Joiner(/*max_dist*/ 1000,
                                /*ratio_to_fragment*/ 10,
                                /*gap_ratio*/ 2));
}

}

