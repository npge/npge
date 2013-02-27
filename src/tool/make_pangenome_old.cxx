/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "process.hpp"
#include "Pipe.hpp"
#include "AddSequences.hpp"
#include "AddBlocks.hpp"
#include "Joiner.hpp"
#include "Connector.hpp"
#include "OverlapsResolver.hpp"
#include "FragmentsExpander.hpp"
#include "BlocksExpander.hpp"
#include "Filter.hpp"
#include "CheckNoOverlaps.hpp"
#include "OutputPipe.hpp"

using namespace bloomrepeats;

class CleanUpOld : public Pipe {
public:
    CleanUpOld() {
        using namespace boost;
        shared_ptr<Filter> filter = make_shared<Filter>();
        filter->set_min_fragment_length(10);
        filter->set_no_options(true);
        add(filter);
        add(new Connector);
        add(filter);
        shared_ptr<OverlapsResolver> resolver = make_shared<OverlapsResolver>();
        add(resolver);
        shared_ptr<Joiner> joiner = make_shared<Joiner>(0);
        joiner->set_no_options(true);
        add(joiner);
        add(filter);
        add(new BlocksExpander);
        add(resolver);
        add(new FragmentsExpander);
        add(new Filter);
        add(new Joiner(/*max_dist*/ 1000,
                                    /*ratio_to_fragment*/ 10,
                                    /*gap_ratio*/ 2));
    }
};

class CleanUpPipe : public Pipe {
public:
    CleanUpPipe() {
        set_empty_block_set();
        add(new AddSequences);
        add(new AddBlocks);
        add(new Connector);
        add(new CleanUpOld);
        add(new CheckNoOverlaps);
        add(new OutputPipe);
    }
};

int main(int argc, char** argv) {
    return process(argc, argv, new CleanUpPipe,
                   "Resolve overlaps and expand blocks");
}

