/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <ostream>

#include "Stats.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "FragmentCollection.hpp"
#include "Rest.hpp"
#include "boundaries.hpp"
#include "block_stat.hpp"
#include "throw_assert.hpp"
#include "report_list.hpp"

namespace npge {

Stats::Stats():
    file_writer_(this, "out-stats", "Output file with statistics") {
    declare_bs("target", "Target blockset");
    add_opt("short-stats", "Print shorter stats", false);
}

// TODO rename Boundaries to smth
typedef Boundaries Integers;

bool fragment_has_overlaps(const VectorFc& fc,
                           Fragment* f) {
    Fragment* next = fc.next(f);
    Fragment* prev = fc.prev(f);
    return (next && next != f && f->common_positions(*next)) ||
           (prev && prev != f && f->common_positions(*prev));
}

template<typename T1, typename T2>
static void report_part(std::ostream& o, const std::string& name,
                        T1 part, T2 total) {
    o << name << ": " << part;
    if (total) {
        Decimal portion = Decimal(part) / Decimal(total);
        Decimal percentage = portion * 100;
        o << " (" << percentage << "%)";
    }
    o << "\n";
}

static void blocks_lengths(std::ostream& out, BlockSetPtr bs) {
    int blocks_sum = 0;
    BOOST_FOREACH (Block* block, *bs) {
        blocks_sum += block->alignment_length();
    }
    out << "Total length of blocks: " << blocks_sum << "\n";
}

static void report_weighted_average_identity(
    std::ostream& out, BlockSetPtr bs) {
    double identity_wsum = 0.0;
    pos_t length_sum = 0;
    BOOST_FOREACH (Block* block, *bs) {
        AlignmentStat al_stat;
        make_stat(al_stat, block);
        double identity = block_identity(al_stat).to_d();
        identity_wsum += double(block->alignment_length()) *
            identity;
        length_sum += block->alignment_length();
    }
    if (!bs->empty()) {
        double identity = identity_wsum / length_sum;
        out << "Identity of joined blocks: ";
        out << identity << "\n";
    }
}

void Stats::run_impl() const {
    int shorter_stats = opt_value("short-stats").as<bool>();
    int blocks_with_alignment = 0, total_fragments = 0;
    int overlap_fragments = 0, overlap_blocks = 0;
    size_t total_nucl = 0, total_seq_length = 0;
    size_t unique_nucl = 0;
    VectorFc fc;
    fc.add_bs(*block_set());
    fc.prepare();
    Integers block_size, fragment_length;
    Decimals spreading; // (max - min) / avg fragment length
    Decimals identity;
    Decimals gc;
    BOOST_FOREACH (Block* b, *block_set()) {
        block_size.push_back(b->size());
        AlignmentStat al_stat;
        make_stat(al_stat, b);
        identity.push_back(block_identity(al_stat));
        gc.push_back(al_stat.gc());
        bool has_overlaps = false;
        BOOST_FOREACH (Fragment* f, *b) {
            total_nucl += f->length();
            fragment_length.push_back(f->length());
            total_fragments += 1;
            if (fragment_has_overlaps(fc, f)) {
                overlap_fragments += 1;
                has_overlaps = true;
            }
        }
        if (has_overlaps) {
            overlap_blocks += 1;
        }
        if (!b->empty() && al_stat.alignment_rows() == b->size()) {
            blocks_with_alignment += 1;
        }
        if (!b->empty()) {
            spreading.push_back(al_stat.spreading());
        }
    }
    Integers seq_length;
    BOOST_FOREACH (SequencePtr s, block_set()->seqs()) {
        seq_length.push_back(s->size());
        total_seq_length += s->size();
    }
    Rest r;
    r.set_other(block_set());
    r.set_empty_block_set();
    r.run();
    Integers unique_length;
    BOOST_FOREACH (Block* b, *r.block_set()) {
        BOOST_FOREACH (Fragment* f, *b) {
            unique_length.push_back(f->length());
            unique_nucl += f->length();
        }
    }
    ASSERT_GTE(total_seq_length, unique_nucl);
    size_t seq_nucl_in_blocks = total_seq_length - unique_nucl;
    int bss = block_set()->size();
    Decimal fpb = bss ? Decimal(total_fragments) / bss : 0;
    std::ostream& out = file_writer_.output();
    if (total_fragments != bss) {
        out << "fragments / blocks = ";
        out << total_fragments << " / " << bss << " = " << fpb << "\n";
        out << "Block identity:";
        report_list(out, identity);
        report_weighted_average_identity(out, block_set());
        if (!shorter_stats) {
            out << "Block sizes:";
            report_list(out, block_size);
        }
    } else {
        out << "fragments = blocks = " << total_fragments << "\n";
    }
    if (!shorter_stats) {
        out << "Fragment lengths:";
        report_list(out, fragment_length);
        if (total_fragments != bss) {
            out << "Fragment length spreading" << "\n";
            out << "  ((max - min) / avg) inside block:";
            report_list(out, spreading);
        }
        out << "GC content:";
        report_list(out, gc);
    }
    report_part(out, "Length of fragments", total_nucl, total_seq_length);
    if (seq_nucl_in_blocks != total_nucl) {
        report_part(out, "Sequence nucleotides in blocks",
                    seq_nucl_in_blocks, total_seq_length);
    }
    if (blocks_with_alignment != 0 && blocks_with_alignment != bss) {
        report_part(out, "Blocks with alignment",
                    blocks_with_alignment, bss);
    }
    if (overlap_fragments != 0) {
        out << "Fragments with overlaps: " << overlap_fragments << "\n";
        out << "Blocks with overlaps: " << overlap_blocks << "\n";
    }
    blocks_lengths(out, block_set());
}

const char* Stats::name_impl() const {
    return "Print human readable summary and statistics about blockset";
}

}

