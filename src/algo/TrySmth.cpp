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

TrySmth::TrySmth() {
    add(new Clear, "target=smth-copy");
    add(new Union, "target=smth-copy other=target");
    add(new MetaProcessor, "prefix|smth-");
    add(new Union, "target=smth-copy other=target");
    add(new Align, "other=target");
    add(new OverlaplessUnion, "target=target other=smth-copy");
    add(new Clear, "target=smth-copy");
}

}

