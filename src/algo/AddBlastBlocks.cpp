/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "AddBlastBlocks.hpp"
#include "SequencesFromOther.hpp"
#include "Clear.hpp"
#include "ConSeq.hpp"
#include "RemoveNames.hpp"
#include "UniqueNames.hpp"
#include "BlastFinder.hpp"
#include "DeConSeq.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"

namespace bloomrepeats {

class FilterSeqs : public Processor {
public:
    FilterSeqs() {
        add_opt("blast-min-length", "min length of blast hit", 100);
        // FIXME 100
        add_opt_rule("blast-min-length >= 0");
    }

protected:
    bool run_impl() const {
        BlockSet& bs = *block_set();
        std::vector<SequencePtr> seqs = bs.seqs();
        int min_length = opt_value("blast-min-length").as<int>();
        BOOST_FOREACH (const SequencePtr& seq, seqs) {
            if (seq->size() < min_length) {
                bs.remove_sequence(seq);
            }
        }
        return true;
    }
};

AddBlastBlocks::AddBlastBlocks(BlockSetPtr source):
    Pipe(source) {
    add(new SequencesFromOther);
    add(new Clear, "target=consensus --clear-seqs:=1");
    add(new ConSeq, "target=consensus other=other");
    add(new RemoveNames, "target=consensus");
    add(new UniqueNames, "target=consensus");
    add(new FilterSeqs, "target=consensus");
    add(new BlastFinder, "target=consensus");
    add(new DeConSeq, "target=target other=consensus");
}

}

