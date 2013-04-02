/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <ostream>
#include <algorithm>

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

static size_t total_length(const BlockSet& bs) {
    size_t result = 0;
    BOOST_FOREACH (Block* b, bs) {
        BOOST_FOREACH (Fragment* f, *b) {
            result += f->length();
        }
    }
    return result;
}

// FIXME buggy
static int fragment_right_overlap(const Fragment* fragment) {
    size_t result = 0;
    const Fragment* f = fragment;
    while (f = f->neighbor(1)) {
        size_t overlap = fragment->common_positions(*f);
        if (!overlap) {
            break;
        }
        result = std::max(result, overlap);
    }
    if (result > fragment->length()) {
        result = fragment->length();
    }
    return result;
}

static int fragment_gc(const Fragment* f) {
    int gc = 0;
    for (int i = 0; i < f->length(); i++) {
        char nucl = f->at(i);
        if (nucl == 'g' || nucl == 'c') {
            gc += 1;
        }
    }
    return gc;
}

template<typename Vector>
static void report_list(std::ostream& o, const Vector& list) {
    o << " size=" << list.size();
    if (!list.empty()) {
        o << " min=" << *std::min_element(list.begin(), list.end());
        o << " avg=" << avg_element(list);
        o << " max=" << *std::max_element(list.begin(), list.end());
    }
    o << std::endl;
}

bool Stats::run_impl() const {
    Connector c;
    c.apply(block_set());
    int blocks_with_alignment = 0, total_fragments = 0;
    int empty_blocks = 0, one_fragment_blocks = 0, one_fragment_nucl = 0;
    int short_fragments = 0, blocks_with_short = 0, small_blocks = 0;
    int overlap_fragments = 0, overlap_blocks = 0;
    size_t overlap_fr_nucl = 0, non_overlap_fr_nucl = 0, overlap_seq_nucl = 0;
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
        if (b->empty()) {
            empty_blocks += 1;
        }
        if (b->size() == 1) {
            one_fragment_blocks += 1;
            one_fragment_nucl += b->front()->length();
        }
        bool has_short = false;
        bool has_alignment = true;
        bool has_overlaps = false;
        BOOST_FOREACH (Fragment* f, *b) {
            total_nucl += f->length();
            fragment_length.push_back(f->length());
            float this_gc = f->length() ? float(fragment_gc(f)) / f->length() : 0;
            gc.push_back(this_gc);
            total_fragments += 1;
            if (!f->row()) {
                has_alignment = false;
            }
            if (f->length() < min_fragment_length()) {
                short_fragments += 1;
                has_short = true;
            }
            int overlaps = fragment_right_overlap(f);
            //overlap_fr_nucl += overlaps * 2;
            if (overlaps) {
                overlap_fragments += 1;
                has_overlaps = true;
            }
        }
        if (b->size() < min_block_size()) {
            small_blocks += 1;
        }
        if (has_short) {
            blocks_with_short += 1;
        }
        if (has_overlaps) {
            overlap_blocks += 1;
        }
        if (!b->empty() && has_alignment) {
            blocks_with_alignment += 1;
        }
        if (!b->empty()) {
            spreading.push_back(al_stat.spreading);
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
    //BOOST_ASSERT(total_nucl >= overlap_fr_nucl);
    //non_overlap_fr_nucl = total_nucl - overlap_fr_nucl;
    BOOST_ASSERT(total_seq_length >= unique_nucl);
    size_t seq_nucl_in_blocks = total_seq_length - unique_nucl;
    //BOOST_ASSERT(seq_nucl_in_blocks >= non_overlap_fr_nucl);
    //overlap_seq_nucl = seq_nucl_in_blocks - non_overlap_fr_nucl;
    output() << "Sequences: " << block_set()->seqs().size() << std::endl;
    output() << "Number of blocks: " << block_set()->size() << std::endl;
    output() << "Number of fragments: " << total_fragments << std::endl;
    size_t bss = block_set()->size();
    float fpb = bss ? float(total_fragments) / bss : 0;
    output() << "Average number of fragments per block: " << fpb << std::endl;
    output() << "Blocks with alignment: " << blocks_with_alignment << std::endl;
    output() << "Empty blocks: " << empty_blocks << std::endl;
    output() << "Blocks of one fragment: " << one_fragment_blocks << std::endl;
    output() << "Nucleotides in blocks of one fragment: "
             << one_fragment_nucl << std::endl;
    output() << "Short (<" << min_fragment_length() << " nucl.) fragments: "
             << short_fragments << std::endl;
    output() << "Blocks with short fragments: "
             << blocks_with_short << std::endl;
    if (min_block_size() > 1 + 1) {
        output() << "Small (<" << min_block_size() << " fragments) blocks: "
                 << small_blocks << std::endl;
    }
    output() << "Block sizes:";
    report_list(output(), block_size);
    output() << "Fragment lengths:";
    report_list(output(), fragment_length);
    output() << "GC content:";
    report_list(output(), gc);
    output() << "Fragment length spreading ((max - min) / avg) inside block:";
    output() << std::endl << "  ";
    report_list(output(), spreading);
    output() << "Block identity:";
    report_list(output(), identity);
    output() << "Sequence lengths:";
    report_list(output(), seq_length);
    output() << "Length of sequences: " << total_seq_length << std::endl;
    output() << "Length of fragments: " << total_nucl << std::endl;
    output() << "Unique fragments lengths:";
    report_list(output(), unique_length);
    output() << "Length of unique fragments: " << unique_nucl << std::endl;
    output() << "Sequence nucleotides in blocks: "
             << seq_nucl_in_blocks << std::endl;
    //output() << "Sequence nucleotides in overlaps: "
    //    << overlap_seq_nucl << std::endl;
    //output() << "Fragment nucleotides in overlaps: "
    //    << overlap_fr_nucl << std::endl;
    //output() << "Non-overlapping nucleotides in blocks: "
    //    << non_overlap_fr_nucl << std::endl;
    output() << "Fragments with overlaps: " << overlap_fragments << std::endl;
    output() << "Blocks with overlaps: " << overlap_blocks << std::endl;
    if (!overlap_fragments && unique_length.empty()) {
        output() << "[pangenome]" << std::endl;
    } else {
        output() << "[not pangenome]" << std::endl;
    }
    return true; // because of Connector
}

const char* Stats::name_impl() const {
    return "Print human readable summary and statistics about block set";
}

}

