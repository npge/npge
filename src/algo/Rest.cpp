/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "Rest.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "FragmentCollection.hpp"
#include "throw_assert.hpp"

namespace npge {

Rest::Rest(const BlockSetPtr& source) {
    set_other(source);
    add_opt("skip-rest",
            "do not add unique fragments to blockset",
            false);
    declare_bs("target", "Where created blocks are added");
    declare_bs("other", "Input blocks");
}

static void add_f(BlockSet& self, Sequence* seq,
                  int min_pos, int max_pos) {
    min_pos = std::max(0, min_pos);
    max_pos = std::min(int(seq->size()) - 1, max_pos);
    if (min_pos > max_pos) {
        return;
    }
    Fragment* new_f = new Fragment(seq, min_pos, max_pos);
    Block* new_b = new Block;
    new_b->insert(new_f);
    self.insert(new_b);
}

void Rest::run_impl() const {
    if (opt_value("skip-rest").as<bool>()) {
        return;
    }
    BlockSet& self = *block_set();
    VectorFc fc;
    fc.add_bs(*other());
    fc.prepare();
    self.add_sequences(other()->seqs());
    std::set<Sequence*> seqs;
    BOOST_FOREACH (SequencePtr s, other()->seqs()) {
        seqs.insert(s.get());
    }
    BOOST_FOREACH (Sequence* s, fc.seqs()) {
        seqs.insert(s);
    }
    BOOST_FOREACH (Sequence* seq, seqs) {
        if (!fc.has_seq(seq)) {
            add_f(self, seq, 0, int(seq->size()) - 1);
        } else {
            const Fragments& ff = fc.fragments_of(seq);
            ASSERT_GT(ff.size(), 0);
            add_f(self, seq, 0, int(ff[0]->min_pos()) - 1);
            for (int i = 1; i < ff.size(); i++) {
                add_f(self, seq,
                      int(ff[i - 1]->max_pos()) + 1,
                      int(ff[i]->min_pos()) - 1);
            }
            add_f(self, seq, int(ff.back()->max_pos()) + 1,
                  int(seq->size()) - 1);
        }
    }
}

const char* Rest::name_impl() const {
    return "Add to target blocks of nucleotides, "
           "not included to other";
}

}

