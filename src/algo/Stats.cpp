/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
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
#include "Connector.hpp"
#include "Rest.hpp"
#include "boundaries.hpp"
#include "block_stat.hpp"
#include "throw_assert.hpp"
#include "report_list.hpp"

namespace npge {

Stats::Stats():
    file_writer_(this, "out-stats", "Output file with statistics") {
    declare_bs("target", "Target blockset");
}

// TODO rename Boundaries to smth
typedef Boundaries Integers;

bool fragment_has_overlaps(const Fragment* f) {
    return (f->next() && f->common_positions(*f->next())) ||
           (f->prev() && f->common_positions(*f->prev()));
}

template<typename T1, typename T2>
static void report_part(std::ostream& o, const std::string& name,
                        T1 part, T2 total) {
    o << name << ": " << part;
    if (total) {
        double portion = float(part) / float(total);
        double percentage = portion * 100;
        o << " (" << percentage << "%)";
    }
    o << "\n";
}

void Stats::run_impl() const {
    Connector c;
    c.apply(block_set());
    int blocks_with_alignment = 0, total_fragments = 0;
    int overlap_fragments = 0, overlap_blocks = 0;
    size_t total_nucl = 0, total_seq_length = 0;
    size_t unique_nucl = 0;
    Integers block_size, fragment_length;
    Floats spreading; // (max - min) / avg fragment length
    Floats identity;
    Floats gc;
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
            if (fragment_has_overlaps(f)) {
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
    size_t bss = block_set()->size();
    double fpb = bss ? float(total_fragments) / bss : 0;
    std::ostream& out = file_writer_.output();
    if (total_fragments != bss) {
        out << "fragments / blocks = ";
        out << total_fragments << " / " << bss << " = " << fpb << "\n";
        out << "Block identity:";
        report_list(out, identity);
        out << "Block sizes:";
        report_list(out, block_size);
    } else {
        out << "fragments = blocks = " << total_fragments << "\n";
    }
    out << "Fragment lengths:";
    report_list(out, fragment_length);
    if (total_fragments != bss) {
        out << "Fragment length spreading" << "\n";
        out << "  ((max - min) / avg) inside block:";
        report_list(out, spreading);
    }
    out << "GC content:";
    report_list(out, gc);
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
}

const char* Stats::name_impl() const {
    return "Print human readable summary and statistics about block set";
}

}

