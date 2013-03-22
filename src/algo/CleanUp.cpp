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

class CleanUpLoop : public Pipe {
public:
    CleanUpLoop() {
        set_max_loops(3);
        FragmentsExpander* expander = new FragmentsExpander;
        expander->set_max_overlap(200);
        add(expander);
        add(new Filter);
        add(new OverlapsResolver2, "target=target other=target");
        add(new Connector);
        add(new Joiner(/*max_dist*/ 100,
                                    /*ratio_to_fragment*/ 0.5,
                                    /*gap_ratio*/ 1.5));
    }
};

CleanUp::CleanUp() {
    add(new CleanUpLoop);
    add(new Filter);
    add(new Connector);
    FragmentsExpander* expander = new FragmentsExpander;
    expander->set_no_options(true);
    add(expander);
}

}

