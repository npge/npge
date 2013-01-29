/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "PrintOverlaps.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

PrintOverlaps::PrintOverlaps():
    print_block_(true),
    print_fragment_(false),
    top_scale_(true),
    bottom_scale_(true),
    width_(76),
    marker_('*') {
    set_prefix("overlaps-");
}

typedef std::vector<const Fragment*> Fragments;

static Fragments overlapping_fragments(const Fragment* fragment) {
    Fragments fragments;
    for (int ori = -1; ori <= 1; ori += 2) {
        const Fragment* f = fragment->neighbor(ori);
        while (f && fragment->common_positions(*f)) {
            fragments.push_back(f);
            f = f->neighbor(ori);
        }
    }
    return fragments;
}

std::string fragment_name(const PrintOverlaps* self, const Fragment* f) {
    std::string result;
    BOOST_ASSERT(f->block());
    if (self->print_block()) {
        result += f->block()->name();
    }
    if (self->print_block() && self->print_fragment()) {
        result += "/";
    }
    if (self->print_fragment()) {
        result += f->id();
    }
    if (!result.empty()) {
        result += " ";
    }
    return result;
}

static void print_scale(const PrintOverlaps* self, std::ostream& o,
                        int name_length, int block_length, Block* block) {
    const int SCALE_STEP = 10; // chars
    o << block->name() << ' ';
    o << std::string(name_length - (block->name().size() + 1), ' ');
    o << ' '; // for '|'
    int diagram_length = self->width() - name_length;
    diagram_length -= 2; // for '|'
    for (int i = 0; i < diagram_length; i += SCALE_STEP) {
        int block_pos = i * block_length / diagram_length;
        std::string block_pos_str = boost::lexical_cast<std::string>(block_pos);
        if (i + block_pos_str.size() <= diagram_length) {
            o << block_pos_str;
            if (i + SCALE_STEP < diagram_length) {
                o << std::string(SCALE_STEP - block_pos_str.size(), ' ');
            }
        }
    }
    o << std::endl;
}

static int block_pos(const Fragment* f, int f_pos, int block_length) {
    if (f->row()) {
        return f->row()->map_to_alignment(f_pos);
    } else {
        return f_pos * block_length / f->length();
    }
}

static void print_overlap(const PrintOverlaps* self, std::ostream& o,
                          const std::string& name, int name_length,
                          int block_length,
                          const Fragment* fragment, const Fragment* f) {
    o << name;
    o << std::string(name_length - name.size(), ' ');
    o << '|';
    Fragment c = fragment->common_fragment(*f);
    int f_begin = std::abs(int(c.begin_pos()) - int(fragment->begin_pos()));
    int delta_last = std::abs(int(c.last_pos()) - int(fragment->last_pos()));
    int f_last = fragment->length() - delta_last - 1;
    int b_begin = block_pos(fragment, f_begin, block_length);
    int b_last = block_pos(fragment, f_last, block_length);
    int diagram_length = self->width() - name_length;
    diagram_length -= 2; // for '|'
    std::string diagram(diagram_length, ' ');
    int d_begin = b_begin * diagram_length / block_length;
    int d_last = b_last * diagram_length / block_length;
    for (int i = d_begin; i <= d_last; i++) {
        diagram[i] = self->marker();
    }
    o << diagram;
    o << '|';
    o << std::endl;
}

// first from block passed to print_block, second from block key of B2Fs
typedef std::pair<const Fragment*, const Fragment*> FragmentPair;
typedef std::vector<FragmentPair> FragmentPairs;
typedef std::map<const Block*, FragmentPairs> B2Fs;
typedef std::map<const Fragment*, std::string> F2Name;

// returns max length of fragment name
static size_t find_fragments_names(const PrintOverlaps* self,
                                   const B2Fs& overlaps, F2Name& f2name) {
    size_t max_name_length = 0;
    BOOST_FOREACH (const B2Fs::value_type& b_and_ff, overlaps) {
        const FragmentPairs& pairs = b_and_ff.second;
        BOOST_FOREACH (const FragmentPair& pair, pairs) {
            const Fragment* f = pair.second;
            f2name[f] = fragment_name(self, f);
            if (f2name[f].size() > max_name_length) {
                max_name_length = f2name[f].size();
            }
        }
    }
    return max_name_length;
}

void PrintOverlaps::print_block(std::ostream& o, Block* block) const {
    B2Fs overlaps;
    BOOST_FOREACH (const Fragment* fragment, *block) {
        BOOST_FOREACH (const Fragment* f, overlapping_fragments(fragment)) {
            overlaps[f->block()].push_back(std::make_pair(fragment, f));
        }
    }
    size_t max_name_length = 0;
    F2Name f2name;
    if (print_block() || print_fragment()) {
        max_name_length = find_fragments_names(this, overlaps, f2name);
    }
    max_name_length = std::max(max_name_length, block->name().size() + 1);
    int diagram_length = width() - max_name_length;
    int block_length = block->alignment_length();
    if (top_scale()) {
        print_scale(this, o, max_name_length, block_length, block);
    }
    bool empty = true;
    BOOST_FOREACH (const B2Fs::value_type& b_and_ff, overlaps) {
        const FragmentPairs& pairs = b_and_ff.second;
        BOOST_FOREACH (const FragmentPair& pair, pairs) {
            const Fragment* fragment = pair.first;
            const Fragment* f = pair.second;
            print_overlap(this, o, f2name[f], max_name_length, block_length,
                          fragment, f);
            empty = false;
        }
    }
    if (empty) {
        o << std::string(max_name_length, ' ');
        o << "No overlaps" << std::endl;
    }
    if (bottom_scale()) {
        print_scale(this, o, max_name_length, block_length, block);
    }
    o << std::endl;
}

void PrintOverlaps::add_options_impl(po::options_description& desc) const {
    AbstractOutput::add_options_impl(desc);
    bloomrepeats::add_unique_options(desc)
    ("print-block", po::value<bool>()->default_value(print_block()),
     "if block name is printed")
    ("print-fragment", po::value<bool>()->default_value(print_fragment()),
     "if fragment name is printed")
    ("top-scale", po::value<bool>()->default_value(top_scale()),
     "if top scale is printed")
    ("bottom-scale", po::value<bool>()->default_value(bottom_scale()),
     "if bottom scale is printed")
    ("width", po::value<int>()->default_value(width()),
     "max allowed line width of output")
    ("marker", po::value<char>()->default_value(marker()),
     "char used to mark fragment")
   ;
}

void PrintOverlaps::apply_options_impl(const po::variables_map& vm) {
    AbstractOutput::apply_options_impl(vm);
    set_print_block(vm["print-block"].as<bool>());
    set_print_fragment(vm["print-fragment"].as<bool>());
    set_top_scale(vm["top-scale"].as<bool>());
    set_bottom_scale(vm["bottom-scale"].as<bool>());
    set_width(vm["width"].as<int>());
    set_marker(vm["marker"].as<char>());
}

const char* PrintOverlaps::name_impl() const {
    return "Print ASCII diagram with all fragments overlapping with a block";
}

}

