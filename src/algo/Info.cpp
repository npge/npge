/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>

#include "Info.hpp"
#include "Filter.hpp"
#include "SizeLimits.hpp"
#include "RemoveNonStem.hpp"
#include "Union.hpp"
#include "Rest.hpp"
#include "Stats.hpp"
#include "Meta.hpp"
#include "FileWriter.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Sequence.hpp"
#include "report_list.hpp"
#include "throw_assert.hpp"

namespace npge {

Info::Info() {
    stats_ = new Stats;
    stats_->set_parent(this);
    declare_bs("target", "Target blockset");
    declare_bs("g-blocks", "g-blocks");
    add_opt("short-stats", "Print shorter stats", false);
}

// TODO rename Boundaries to smth
typedef Boundaries Integers;

static int blocks_lengths(std::ostream& out, BlockSetPtr bs) {
    int minor_sum = 0;
    int regular_sum = 0;
    BOOST_FOREACH (Block* block, *bs) {
        ASSERT_GT(block->name().length(), 0);
        if (block->name()[0] == 'm') {
            minor_sum += block->alignment_length();
        } else {
            regular_sum += block->alignment_length();
        }
    }
    out << "Blocks' lengths: regular + minor = ";
    out << regular_sum << " + " << minor_sum << " = ";
    out << (regular_sum + minor_sum) << "\n";
    return regular_sum + minor_sum;
}

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
    //
    int npg_length = blocks_lengths(out, block_set());
    stats_->set_npg_length(npg_length);
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
    return u.block_set();
}

void Info::print_all() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nAll regular blocks of at least 2 fragments:\n";
    BlockSetPtr bs = filter_blocks();
    meta()->get("RemoveMinorBlocks")->apply(bs);
    stats_->apply(bs);
}

static BlockSetPtr filter_by_letter(
    BlockSetPtr block_set,
    char letter) {
    BlockSetPtr bs = new_bs();
    bs->add_sequences(block_set->seqs());
    BOOST_FOREACH (Block* block, *block_set) {
        std::string name = block->name();
        if (!name.empty() && name[0] == letter) {
            bs->insert(block->clone());
        }
    }
    return bs;
}

void Info::print_rest() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nRest (blocks of 1 fragment but not minor):\n";
    BlockSetPtr bs = filter_by_letter(block_set(), 'u');
    stats_->apply(bs);
}

void Info::print_minor() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nMinor blocks (too short to say smth about):\n";
    BlockSetPtr bs = filter_by_letter(block_set(), 'm');
    stats_->apply(bs);
}

void Info::print_hemi() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nPartial blocks ";
    out << "(represented once in subset of genomes):\n";
    BlockSetPtr bs = filter_by_letter(block_set(), 'h');
    stats_->apply(bs);
}

void Info::print_repeats() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nBlocks with repeats ";
    out << "(at least two copies in at least one genome):\n";
    BlockSetPtr bs = filter_by_letter(block_set(), 'r');
    stats_->apply(bs);
}

void Info::print_stem() const {
    std::ostream& out = stats_->file_writer().output();
    out << "\nExact stem blocks (represented in all genomes) "
        "but not minor:\n";
    BlockSetPtr bs = filter_blocks();
    meta()->get("RemoveMinorBlocks")->apply(bs);
    RemoveNonStem stem;
    stem.set_opt_value("exact", true);
    stem.set_block_set(bs);
    try {
        stem.run();
        stats_->apply(stem.block_set());
    } catch (...) {
        out << "\nWarning: failed to build stem\n";
    }
    out << "\n";
}

void Info::print_global() const {
    BlockSetPtr g_bs = get_bs("g-blocks");
    if (g_bs->empty()) {
        return;
    }
    std::ostream& out = stats_->file_writer().output();
    out << "\nG-blocks:\n";
    stats_->apply(g_bs);
    out << "\n";
}

void Info::run_impl() const {
    int shorter_stats = opt_value("short-stats").as<bool>();
    print_seq();
    if (!shorter_stats) {
        print_all();
    }
    print_stem();
    if (!shorter_stats) {
        print_hemi();
        print_repeats();
        print_rest();
        print_minor();
        print_global();
    }
}

const char* Info::name_impl() const {
    return "Print human readable summary and statistics";
}

}

