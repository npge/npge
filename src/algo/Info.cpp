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
    add_opt("short-stats", "Print shorter stats", false);
}

// TODO rename Boundaries to smth
typedef Boundaries Integers;

void Info::print_seq() const {
    std::ostream& out = stats_->file_writer().output();
    pos_t total_seq_length = 0;
    Integers seq_length;
    typedef std::map<std::string, int> Genome2Length;
    Genome2Length g2l;
    BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
        seq_length.push_back(s->size());
        g2l[s->genome()] += s->size();
        total_seq_length += s->size();
    }
    out << "Number of sequences: " << block_set()->seqs().size() << "\n";
    out << "Sequence lengths:";
    report_list(out, seq_length);
    out << "Total length of sequences: " << total_seq_length << std::endl;
    Integers genomes_length;
    BOOST_FOREACH (const Genome2Length::value_type& kv, g2l) {
        genomes_length.push_back(kv.second);
    }
    out << "Genomes:";
    report_list(out, genomes_length);
}

BlockSetPtr Info::filter_blocks() const {
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
    return u.block_set();
}

void Info::print_all() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nAll non-minor blocks of at least 2 fragments:\n";
    BlockSetPtr bs = filter_blocks();
    stats_->apply(bs);
}

void Info::print_rest() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nRest (sequence parts not covered by blocks of >= 2 fr.):\n";
    BlockSetPtr bs = filter_blocks();
    Rest rest;
    rest.set_other(bs);
    rest.run();
    rest.block_set()->add_sequences(block_set()->seqs());
    stats_->apply(rest.block_set());
}

void Info::print_stem() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nExact stem blocks (represented in all genomes) "
        "but not minor:\n";
    BlockSetPtr bs = filter_blocks();
    Stem stem;
    stem.set_opt_value("exact", true);
    stem.set_block_set(bs);
    try {
        stem.run();
        stats_->apply(stem.block_set());
    } catch (...) {
        out << "\nFailed to build stem\n";
    }
    out << "\n";
}

void Info::run_impl() const {
    int shorter_stats = opt_value("short-stats").as<bool>();
    if (!shorter_stats) {
        print_seq();
        print_all();
        print_rest();
    }
    print_stem();
}

const char* Info::name_impl() const {
    return "Print human readable summary and statistics";
}

}

