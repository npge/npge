/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "CleanUp.hpp"
#include "Joiner.hpp"
#include "Connector.hpp"
#include "OverlapsResolver2.hpp"
#include "FragmentsExpander.hpp"
#include "Filter.hpp"

namespace bloomrepeats {

CleanUp::CleanUp() {
    using namespace boost;
    FragmentsExpander* expander = new FragmentsExpander;
    expander->set_max_overlap(200);
    add(expander);
    add(new Filter);
    add(new OverlapsResolver2, "target=target other=target");
    add(new Connector);
    add(new Joiner(/*max_dist*/ 1000,
                                /*ratio_to_fragment*/ 10,
                                /*gap_ratio*/ 2));
}

}

