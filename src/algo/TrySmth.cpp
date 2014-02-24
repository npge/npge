/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "TrySmth.hpp"
#include "Union.hpp"
#include "MetaProcessor.hpp"
#include "Clear.hpp"
#include "OverlaplessUnion.hpp"
#include "Align.hpp"

namespace bloomrepeats {

class AddingLoop : public Pipe {
public:
    AddingLoop() {
        set_max_loops(-1);
        add(new OverlaplessUnion);
        add(new Align);
    }
};

TrySmth::TrySmth() {
    add(new Union, "target=smth-copy other=target");
    add(new MetaProcessor, "prefix|smth-");
    add(new Union, "target=smth-copy other=target");
    add(new Clear, "target=target");
    add(new AddingLoop, "target=target other=smth-copy");
    add(new Clear, "target=smth-copy --clear-seqs=1");
}

}

