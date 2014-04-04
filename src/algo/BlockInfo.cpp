/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <map>
#include <boost/foreach.hpp>

#include "BlockInfo.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "block_stat.hpp"
#include "global.hpp"

namespace bloomrepeats {

BlockInfo::BlockInfo(const std::string& prefix) {
    set_opt_prefix(prefix);
    add_opt("count-seqs",
            "Add columns with orrurences in each sequence",
            false);
    declare_bs("target", "Target blockset");
}

struct CmpSeqs {
    bool operator()(const SequencePtr& a, const SequencePtr& b) const {
        return a->name() < b->name();
    }
};

void BlockInfo::print_header(std::ostream& o) const {
    o << "block";
    o << '\t' << "fragments";
    o << '\t' << "cols";
    o << '\t' << "ident-nogap";
    o << '\t' << "ident-gap";
    o << '\t' << "noident-nogap";
    o << '\t' << "noident-gap";
    o << '\t' << "pure-gap";
    o << '\t' << "ident";
    o << '\t' << "GC";
    bool count_seqs = opt_value("count-seqs").as<bool>();
    if (count_seqs) {
        std::vector<SequencePtr> seqs = block_set()->seqs();
        std::sort(seqs.begin(), seqs.end(), CmpSeqs());
        BOOST_FOREACH (const SequencePtr& seq, seqs) {
            o << '\t' << seq->name();
        }
    }
    o << "\n";
}

void BlockInfo::print_block(std::ostream& o, Block* block) const {
    o << block->name();
    o << '\t' << block->size();
    AlignmentStat stat;
    make_stat(stat, block);
    o << '\t' << stat.total();
    o << '\t' << stat.ident_nogap();
    o << '\t' << stat.ident_gap();
    o << '\t' << stat.noident_nogap();
    o << '\t' << stat.noident_gap();
    o << '\t' << stat.pure_gap();
    o << '\t' << block_identity(stat);
    o << '\t' << stat.gc();
    bool count_seqs = opt_value("count-seqs").as<bool>();
    if (count_seqs) {
        std::vector<SequencePtr> seqs = block_set()->seqs();
        std::sort(seqs.begin(), seqs.end(), CmpSeqs());
        std::map<Sequence*, int> occurences;
        BOOST_FOREACH (Fragment* f, *block) {
            occurences[f->seq()] += 1;
        }
        BOOST_FOREACH (const SequencePtr& seq, seqs) {
            o << '\t' << occurences[seq.get()];
        }
    }
    o << "\n";
}

const char* BlockInfo::name_impl() const {
    return "Print information about blocks";
}

}

