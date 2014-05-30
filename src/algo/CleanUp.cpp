/*
 * NPG-explorer, Nucleotide PanGenome explorer
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
#include "to_s.hpp"

namespace npge {

class CleanUpLoop : public Pipe {
public:
    CleanUpLoop() {
        set_max_loops(3);
        add(new FragmentsExpander,
            "--max-overlap=$EXPANDER_MAX_OVERLAP");
        add(new Filter);
        add(new OverlapsResolver2, "target=target other=target");
        add(new Connector);
        add(new Joiner, "--join-max-dist=$JOINER_MAX_DIST "
            "--join-to-fragment=$JOINER_RATIO_TO_FRAGMENT "
            "--join-to-gap=$JOINER_GAP_RATIO");
    }
};

CleanUp::CleanUp() {
    add(new CleanUpLoop);
    add(new Filter);
    add(new Connector);
    add(new FragmentsExpander, "no_options");
    declare_bs("target", "Target blockset");
}

const char* CleanUp::name_impl() const {
    return "Connect, resolve overlaps, expand, filter (deprecated)";
}

}

