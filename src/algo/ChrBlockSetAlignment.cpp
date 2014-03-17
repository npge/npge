/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <vector>
#include <boost/foreach.hpp>

#include "ChrBlockSetAlignment.hpp"
#include "BlockSetAlignment.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"

namespace bloomrepeats {

typedef std::vector<std::string> SeqGroups;

ChrBlockSetAlignment::ChrBlockSetAlignment() {
    bsa_ = new BlockSetAlignment;
    bsa_->set_parent(this);
    bsa_->point_bs("target=target", this);
    BOOST_ASSERT(bsa_->block_set() == block_set());
}

bool ChrBlockSetAlignment::run_impl() const {
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
    return false;
}

const char* ChrBlockSetAlignment::name_impl() const {
    return "Build block set alignments from all chromosomes";
}

}

