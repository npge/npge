/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <ostream>
#include <algorithm>
#include <boost/format.hpp>

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

namespace bloomrepeats {

void Stats::add_options_impl(po::options_description& desc) const {
    SizeLimits::add_options_impl(desc);
    add_unique_options(desc)
    ("out-stats", po::value<std::string>(), "Output file with statistics")
   ;
}

void Stats::apply_options_impl(const po::variables_map& vm) {
    SizeLimits::apply_options_impl(vm);
    if (vm.count("out-stats")) {
        set_output_file(vm["out-stats"].as<std::string>());
    }
}

// TODO rename Boundaries to smth
typedef Boundaries Integers;

bool fragment_has_overlaps(const Fragment* f) {
    return (f->next() && f->common_positions(*f->next())) ||
           (f->prev() && f->common_positions(*f->prev()));
}

template<typename Vector>
static void report_list(std::ostream& o, const Vector& list) {
    o << " number=" << list.size();
    if (!list.empty()) {
        typename Vector::value_type min_value, max_value;
        min_value = *std::min_element(list.begin(), list.end());
        max_value = *std::max_element(list.begin(), list.end());
        if (min_value == max_value) {
            o << " all=" << min_value;
        } else {
            o << " min=" << min_value;
            o << " median=" << median_element(list);
            double avg = avg_element_double(list);
            boost::format double_2("%.2f");
            o << " avg=" << str(double_2 % avg);
            o << " max=" << max_value;
        }
    }
    o << std::endl;
}

template<typename T1, typename T2>
static void report_part(std::ostream& o, const std::string& name,
                        T1 part, T2 total) {
    o << name << ": " << part;
    float portion = float(part) / float(total);
    float percentage = portion * 100;
    o << " (" << percentage << "%)\n";
}

bool Stats::run_impl() const {
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
    BOOST_ASSERT(total_seq_length >= unique_nucl);
    size_t seq_nucl_in_blocks = total_seq_length - unique_nucl;
    size_t bss = block_set()->size();
    float fpb = bss ? float(total_fragments) / bss : 0;
    output() << "fragments / blocks = ";
    output() << total_fragments << " / " << bss << " = " << fpb << "\n";
    output() << "Block identity:";
    report_list(output(), identity);
    output() << "Block sizes:";
    report_list(output(), block_size);
    output() << "Fragment lengths:";
    report_list(output(), fragment_length);
    output() << "Fragment length spreading ((max - min) / avg) inside block:";
    output() << std::endl << "  ";
    report_list(output(), spreading);
    output() << "GC content:";
    report_list(output(), gc);
    report_part(output(), "Length of fragments", total_nucl, total_seq_length);
    if (seq_nucl_in_blocks != total_nucl) {
        report_part(output(), "Sequence nucleotides in blocks",
                    seq_nucl_in_blocks, total_seq_length);
    }
    if (blocks_with_alignment != 0 && blocks_with_alignment != bss) {
        report_part(output(), "Blocks with alignment",
                    blocks_with_alignment, bss);
    }
    if (overlap_fragments != 0) {
        output() << "Fragments with overlaps: " << overlap_fragments << "\n";
        output() << "Blocks with overlaps: " << overlap_blocks << "\n";
    }
    return false;
}

const char* Stats::name_impl() const {
    return "Print human readable summary and statistics about block set";
}

}

