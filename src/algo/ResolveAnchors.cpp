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
#include "OverlapsResolver2.hpp"
#include "ConSeq.hpp"
#include "DeConSeq.hpp"

namespace bloomrepeats {

ResolveAnchors::ResolveAnchors(BlockSetPtr source):
    Pipe(source) {
    add(new Connector, "target=other");
    add(new Rest, "target=other other=other");
    set_bs("cons", new_bs());
    add(new ConSeq, "target=cons other=other");
    set_bs("anchors", new_bs());
    add(new AnchorFinder, "target=anchors other=cons");
    add(new OverlapsResolver2, "target=anchors other=anchors");
    add(new Connector, "target=anchors");
    add(new Rest, "target=anchors other=anchors");
    add(new SequencesFromOther, "target=target other=other");
    add(new DeConSeq, "target=target other=anchors");
}

}

