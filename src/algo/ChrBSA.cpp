/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <vector>
#include <boost/foreach.hpp>

#include "ChrBSA.hpp"
#include "FindBSA.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace bloomrepeats {

typedef Strings SeqGroups;

ChrBSA::ChrBSA() {
    bsa_ = new FindBSA;
    bsa_->set_parent(this);
    bsa_->point_bs("target=target", this);
    ASSERT_EQ(bsa_->block_set(), block_set());
    declare_bs("target", "Target blockset");
}

void ChrBSA::run_impl() const {
    BlockSet& bs = *block_set();
    std::set<std::string> chrs;
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        chrs.insert(seq->chromosome());
    }
    BOOST_FOREACH (std::string chr, chrs) {
        SeqGroups seq_groups;
        seq_groups.push_back(chr);
        bsa_->set_opt_value("bsa-seqs", seq_groups);
        bsa_->set_opt_value("bsa-name", chr);
        bsa_->run();
    }
}

const char* ChrBSA::name_impl() const {
    return "Build block set alignments from all chromosomes";
}

}

