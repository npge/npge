/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "AddBlastBlocks.hpp"
#include "SequencesFromOther.hpp"
#include "Clear.hpp"
#include "Union.hpp"
#include "SliceNless.hpp"
#include "ConSeq.hpp"
#include "BlastFinder.hpp"
#include "DeConSeq.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "rand_name.hpp"

namespace npge {

class FilterSeqs : public Processor {
public:
    FilterSeqs() {
        add_gopt("blast-min-length", "min length of blast hit",
                 "MIN_LENGTH");
        add_opt_rule("blast-min-length >= 0");
    }

protected:
    void run_impl() const {
        BlockSet& bs = *block_set();
        std::vector<SequencePtr> seqs = bs.seqs();
        int min_length = opt_value("blast-min-length").as<int>();
        std::set<std::string> names;
        BOOST_FOREACH (const SequencePtr& seq, seqs) {
            seq->set_description("");
            if (seq->size() < min_length) {
                bs.remove_sequence(seq);
            } else {
                if (seq->name().length() > 16) {
                    seq->set_name("");
                }
                while (seq->name().empty() ||
                        names.find(seq->name()) != names.end()) {
                    const int RAND_SEQ_NAME_LENGTH = 8;
                    seq->set_name(rand_name(RAND_SEQ_NAME_LENGTH));
                }
                names.insert(seq->name());
            }
        }
    }
};

AddBlastBlocks::AddBlastBlocks() {
    add(new SequencesFromOther);
    add(new Union, "target=slicenless other=other");
    add(new SliceNless, "target=slicenless");
    add(new ConSeq, "target=consensus other=slicenless");
    add(new FilterSeqs, "target=consensus");
    add(new BlastFinder, "target=consensus");
    add(new DeConSeq, "target=target other=consensus");
    add(new Clear, "target=slicenless --clear-seqs:=1");
    add(new Clear, "target=consensus --clear-seqs:=1");
    declare_bs("other", "Input blocks");
    declare_bs("target", "Blast hits sliced from input blocks");
}

const char* AddBlastBlocks::name_impl() const {
    return "Add blast hits as blocks";
}

}

