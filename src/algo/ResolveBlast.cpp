/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "ResolveBlast.hpp"
#include "SequencesFromOther.hpp"
#include "Rest.hpp"
#include "AddBlastBlocks.hpp"
#include "OverlapsResolver2.hpp"
#include "ConSeq.hpp"
#include "DeConSeq.hpp"

namespace bloomrepeats {

ResolveBlast::ResolveBlast(BlockSetPtr source):
    Pipe(source) {
    add(new Rest, "target=other other=other");
    set_bs("cons", new_bs());
    add(new ConSeq, "target=cons other=other");
    add(new Rest, "target=cons other=cons");
    add(new AddBlastBlocks, "target=cons other=cons");
    add(new OverlapsResolver2, "target=cons other=cons");
    add(new SequencesFromOther, "target=target other=other");
    add(new DeConSeq, "target=target other=cons");
}

}

