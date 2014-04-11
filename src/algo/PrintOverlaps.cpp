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
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "PrintOverlaps.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "proportion.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"
#include "to_s.hpp"
#include "global.hpp"

namespace bloomrepeats {

static bool check_marker_length(const Processor* p,
                                std::string& message) {
    if (p->opt_value("marker").as<std::string>().size() != 1) {
        message = "marker must be of one char length";
        return false;
    }
    return true;
}

PrintOverlaps::PrintOverlaps() {
    set_opt_prefix("overlaps-");
    add_opt("print-block", "if block name is printed", true);
    add_opt("print-fragment", "if fragment name is printed", false);
    add_opt("top-scale", "if top scale is printed", true);
    add_opt("bottom-scale", "if bottom scale is printed", true);
    add_opt("width", "max allowed line width of output", 76);
    add_opt("marker", "char used to mark fragment", std::string("*"));
    add_opt_rule("width > 0");
    add_opt_check(boost::bind(check_marker_length, this, _1));
    declare_bs("target", "Target blockset");
}

static Fragments overlapping_fragments(const Fragment* fragment) {
    Fragments fragments;
    for (int ori = -1; ori <= 1; ori += 2) {
        Fragment* f = fragment->neighbor(ori);
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
    bool p_block = self->opt_value("print-block").as<bool>();
    bool p_fragment = self->opt_value("print-fragment").as<bool>();
    if (p_block) {
        result += f->block()->name();
    }
    if (p_block && p_fragment) {
        result += "/";
    }
    if (p_fragment) {
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
    int width = self->opt_value("width").as<int>();
    int diagram_length = width - name_length;
    diagram_length -= 2; // for '|'
    for (int i = 0; i < diagram_length; i += SCALE_STEP) {
        int block_pos = proportion(i, diagram_length, block_length);
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

static void print_overlap(const PrintOverlaps* self, std::ostream& o,
                          const std::string& name, int name_length,
                          int block_length,
                          const Fragment* fragment, const Fragment* f) {
    o << name;
    BOOST_ASSERT(name_length - name.size() >= 0);
    o << std::string(name_length - name.size(), ' ');
    o << '|';
    Fragment c = fragment->common_fragment(*f);
    int f_begin = std::abs(int(c.begin_pos()) - int(fragment->begin_pos()));
    int delta_last = std::abs(int(c.last_pos()) - int(fragment->last_pos()));
    int f_last = fragment->length() - delta_last - 1;
    int b_begin = block_pos(fragment, f_begin, block_length);
    int b_last = block_pos(fragment, f_last, block_length);
    BOOST_ASSERT(0 <= b_begin && b_begin < block_length);
    BOOST_ASSERT(0 <= b_last && b_last < block_length);
    int width = self->opt_value("width").as<int>();
    char marker = self->opt_value("marker").as<std::string>()[0];
    int diagram_length = width - name_length;
    diagram_length -= 2; // for '|'
    BOOST_ASSERT(diagram_length >= 0);
    std::string diagram(diagram_length, ' ');
    int d_begin = proportion(b_begin, block_length, diagram_length);
    int d_last = proportion(b_last, block_length, diagram_length);
    BOOST_ASSERT(0 <= d_begin && d_begin < diagram_length);
    BOOST_ASSERT_MSG(0 <= d_last && d_last < diagram_length,
                     (TO_S(d_last) + " " +
                      TO_S(diagram_length)).c_str());
    for (int i = d_begin; i <= d_last; i++) {
        diagram[i] = marker;
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
    bool p_block = opt_value("print-block").as<bool>();
    bool p_fragment = opt_value("print-fragment").as<bool>();
    bool top_scale = opt_value("top-scale").as<bool>();
    bool bottom_scale = opt_value("bottom-scale").as<bool>();
    int width = opt_value("width").as<int>();
    if (p_block || p_fragment) {
        max_name_length = find_fragments_names(this, overlaps, f2name);
    }
    max_name_length = std::max(max_name_length, block->name().size() + 1);
    int diagram_length = width - max_name_length;
    int block_length = block->alignment_length();
    if (top_scale) {
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
    if (bottom_scale) {
        print_scale(this, o, max_name_length, block_length, block);
    }
    o << std::endl;
}

const char* PrintOverlaps::name_impl() const {
    return "Print ASCII diagram with all fragments overlapping with a block";
}

}

