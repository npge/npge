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
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "boundaries.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

void Stats::add_options_impl(po::options_description& desc) const {
    SizeLimits::add_options_impl(desc);
}

void Stats::apply_options_impl(const po::variables_map& vm) {
    SizeLimits::apply_options_impl(vm);
}

// TODO rename Boundaries to smth
typedef Boundaries Integers;

typedef std::vector<float> Floats;

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

float avg_element(const Floats& floats) {
    BOOST_ASSERT(!floats.empty());
    float sum = 0;
    BOOST_FOREACH (float f, floats) {
        sum += f;
    }
    return sum / floats.size();
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
    int blocks_with_alignment = 0, total_fragments = 0;
    int empty_blocks = 0, one_fragment_blocks = 0;
    int short_fragments = 0, blocks_with_short = 0, small_blocks = 0;
    Integers block_size, fragment_length;
    Floats spreading; // (max - min) / avg fragment length
    Floats identity;
    Floats gc;
    BOOST_FOREACH (Block* b, *block_set()) {
        block_size.push_back(b->size());
        identity.push_back(b->identity());
        if (b->empty()) {
            empty_blocks += 1;
        }
        if (b->size() == 1) {
            one_fragment_blocks += 1;
        }
        Integers lengths;
        bool has_short = false;
        bool has_alignment = true;
        BOOST_FOREACH (Fragment* f, *b) {
            lengths.push_back(f->length());
            fragment_length.push_back(f->length());
            gc.push_back(float(fragment_gc(f)) / f->length());
            total_fragments += 1;
            if (!f->row()) {
                has_alignment = false;
            }
            if (f->length() < min_fragment_length()) {
                short_fragments += 1;
                has_short = true;
            }
        }
        if (b->size() < min_block_size()) {
            small_blocks += 1;
        }
        if (has_short) {
            blocks_with_short += 1;
        }
        if (!b->empty() && has_alignment) {
            blocks_with_alignment += 1;
        }
        if (!lengths.empty()) {
            int max_length = *std::max_element(lengths.begin(), lengths.end());
            int min_length = *std::min_element(lengths.begin(), lengths.end());
            int avg_length = avg_element(lengths);
            spreading.push_back(float(max_length - min_length) / avg_length);
        }
    }
    output() << "Number of blocks: " << block_set()->size() << std::endl;
    output() << "Number of fragments: " << total_fragments << std::endl;
    float fpb = float(total_fragments) / block_set()->size();
    output() << "Average number of fragments per block: " << fpb << std::endl;
    output() << "Blocks with alignment: " << blocks_with_alignment << std::endl;
    output() << "Empty blocks: " << empty_blocks << std::endl;
    output() << "Blocks of one fragment: " << one_fragment_blocks << std::endl;
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
    return false;
}

const char* Stats::name_impl() const {
    return "Print human readable summary and statistics about block set";
}

}

