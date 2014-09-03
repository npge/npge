/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>

#include "Info.hpp"
#include "Filter.hpp"
#include "SizeLimits.hpp"
#include "Stem.hpp"
#include "Union.hpp"
#include "Rest.hpp"
#include "Stats.hpp"
#include "Meta.hpp"
#include "FileWriter.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "report_list.hpp"

namespace npge {

Info::Info() {
    stats_ = new Stats;
    stats_->set_parent(this);
    declare_bs("target", "Target blockset");
}

// TODO rename Boundaries to smth
typedef Boundaries Integers;

void Info::run_impl() const {
    size_t total_seq_length = 0;
    Integers seq_length;
    BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
        seq_length.push_back(s->size());
        total_seq_length += s->size();
    }
    std::ostream& out = stats_->file_writer().output();
    out << "Number of sequences: " << block_set()->seqs().size() << "\n";
    out << "Sequence lengths:";
    report_list(out, seq_length);
    out << "Total length of sequences: " << total_seq_length << std::endl;
    //
    out << "\nAll non-minor blocks of at least 2 fragments:\n";
    Union u;
    u.set_other(block_set());
    u.run();
    u.block_set()->add_sequences(block_set()->seqs());
    LiteFilter filter;
    filter.set_block_set(u.block_set());
    filter.set_opt_value("min-block", 2);
    filter.set_opt_value("min-fragment", 0);
    filter.run();
    meta()->get("RemoveMinorBlocks")->apply(u.block_set());
    stats_->apply(filter.block_set());
    //
    out << "\nRest (sequence parts not covered by blocks of >= 2 fr.):\n";
    Rest rest;
    rest.set_other(filter.block_set());
    rest.run();
    rest.block_set()->add_sequences(block_set()->seqs());
    stats_->apply(rest.block_set());
    //
    out << "\nStem (blocks represented in all genomes) "
        "but not minor:\n";
    Stem stem;
    stem.set_block_set(u.block_set()); // reuse
    try {
        stem.run();
        stats_->apply(stem.block_set());
    } catch (...) {
        out << "\nFailed to build stem\n";
    }
    out << "\n";
}

const char* Info::name_impl() const {
    return "Print human readable summary and statistics";
}

}

