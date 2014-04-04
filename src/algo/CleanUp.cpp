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
#include "config.hpp"
#include "to_s.hpp"

namespace bloomrepeats {

class CleanUpLoop : public Pipe {
public:
    CleanUpLoop() {
        set_max_loops(3);
        add(new FragmentsExpander,
            "--max-overlap=" + TO_S(EXPANDER_MAX_OVERLAP));
        add(new Filter);
        add(new OverlapsResolver2, "target=target other=target");
        add(new Connector);
        add(new Joiner(JOINER_MAX_DIST, JOINER_RATIO_TO_FRAGMENT,
                       JOINER_GAP_RATIO));
    }
};

CleanUp::CleanUp() {
    add(new CleanUpLoop);
    add(new Filter);
    add(new Connector);
    add(new FragmentsExpander, "no_options");
    declare_bs("target", "Target blockset");
}

}

