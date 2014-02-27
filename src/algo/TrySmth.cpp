/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "TrySmth.hpp"
#include "Union.hpp"
#include "MetaProcessor.hpp"
#include "Clear.hpp"
#include "OverlaplessUnion.hpp"
#include "Align.hpp"
#include "Filter.hpp"
#include "UniqueNames.hpp"
#include "RemoveNames.hpp"
#include "SizeLimits.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"

namespace bloomrepeats {

class AddingLoop : public Pipe {
public:
    AddingLoop() {
        set_max_loops(-1);
        add(new OverlaplessUnion);
        add(new Align);
    }
};

class AddingLoopBySize : public Processor {
public:
    AddingLoopBySize() {
        al_ = new AddingLoop;
        al_->set_parent(this);
    }

protected:
    bool run_impl() const {
        std::set<int> sizes;
        BOOST_FOREACH (const Block* block, *other()) {
            sizes.insert(block->size());
        }
        Filter fil;
        allow_everything(&fil);
        fil.set_opt_value("find-subblocks", false);
        fil.set_opt_value("good-to-other", true);
        fil.set_bs("target", other());
        al_->set_bs("target", block_set());
        al_->set_bs("other", fil.other());
        BOOST_REVERSE_FOREACH (int size, sizes) {
            // set is ordered
            fil.other()->clear();
            fil.set_opt_value("min-block", size);
            fil.set_opt_value("max-block", size);
            fil.run();
            al_->run();
        }
    }

private:
    AddingLoop* al_;
};

TrySmth::TrySmth() {
    add(new Union, "target=smth-copy other=target");
    add(new UniqueNames, "target=smth-copy");
    add(new MetaProcessor, "prefix|smth-");
    add(new RemoveNames, "target=target --remove-seqs-names:=0 "
        " --remove-blocks-names:=1");
    add(new Filter, "target=target");
    add(new Union, "target=smth-copy other=target");
    add(new Clear, "target=target");
    add(new AddingLoopBySize, "target=target other=smth-copy");
    add(new UniqueNames, "target=target");
    add(new Clear, "target=smth-copy --clear-seqs:=1");
}

}

