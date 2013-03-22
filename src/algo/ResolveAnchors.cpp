/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "ResolveAnchors.hpp"
#include "SequencesFromOther.hpp"
#include "Connector.hpp"
#include "Rest.hpp"
#include "AnchorFinder.hpp"
#include "CleanUp.hpp"
#include "ConSeq.hpp"
#include "DeConSeq.hpp"

namespace bloomrepeats {

ResolveAnchors::ResolveAnchors(BlockSetPtr source):
    Pipe(source) {
    add(new Connector, "target=other");
    add(new Rest, "target=other other=other");
    set_bs("cons", new_bs());
    add(new ConSeq, "target=cons other=other");
    set_bs("", new_bs());
    add(new AnchorFinder, "target=cons");
    add(new CleanUp, "target=cons other=cons");
    add(new Connector, "target=cons");
    add(new Rest, "target=cons other=cons");
    add(new SequencesFromOther, "target=target other=other");
    add(new DeConSeq, "target=target other=cons");
}

}

