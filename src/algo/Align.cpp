/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "Align.hpp"
#include "ExternalAligner.hpp"
#include "MoveGaps.hpp"
#include "CutGaps.hpp"

namespace bloomrepeats {

Align::Align() {
    add(new ExternalAligner);
    add(new MoveGaps);
    add(new CutGaps);
}

bool Align::apply_to_block(Block* block) {
    BOOST_FOREACH (Processor* p, processors()) {
        BlocksJobs* jobs = boost::polymorphic_downcast<BlocksJobs*>(p);
        jobs->apply_to_block(block);
    }
}

}

