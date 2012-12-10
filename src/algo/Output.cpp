/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "Output.hpp"
#include "Exception.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

Output::Output():
    export_alignment_(true)
{ }

void Output::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("out-file,o", po::value<std::string>()->default_value(file()),
     "output file with all blocks")
    ("out-mask", po::value<std::string>()->default_value(mask()),
     "mask of output files (${block} is replaced with block name)")
    ("export-alignment", po::value<bool>()->default_value(export_alignment()),
     "use alignment information if available")
   ;
}

void Output::apply_options_impl(const po::variables_map& vm) {
    set_file(vm["out-file"].as<std::string>());
    set_mask(vm["out-mask"].as<std::string>());
    if (!file().empty() && !mask().empty()) {
        throw Exception("both 'out-file' and 'out-mask' were specified");
    }
    if (!mask().empty() && mask().find("${block}") == std::string::npos) {
        throw Exception("'out-mask' must contain '${block}'");
    }
    set_export_alignment(vm["export-alignment"].as<bool>());
}

static struct BlockCompareName2 {
    bool operator()(const Block* b1, const Block* b2) const {
        typedef boost::tuple<int, const std::string&> Tie;
        return Tie(-b1->size(), b1->name()) < Tie(-b2->size(), b2->name());
    }
} bcn2;

bool Output::run_impl() const {
    std::vector<Block*> blocks(block_set()->begin(), block_set()->end());
    if (mask().empty()) {
        std::sort(blocks.begin(), blocks.end(), bcn2);
    }
    std::ostream* out = 0;
    if (!file().empty()) {
        out = new std::ofstream(file().c_str());
    } else if (file().empty() && mask().empty()) {
        out = &std::cout;
    }
    BOOST_FOREACH (Block* b, *block_set()) {
        std::ostream* o = out;
        if (!out) {
            using namespace boost::algorithm;
            std::string path = replace_all_copy(mask(), "${block}", b->name());
            o = new std::ofstream(path.c_str());
        }
        BOOST_ASSERT(o);
        BOOST_FOREACH (Fragment* fr, *b) {
            (*o) << '>';
            fr->print_header(*o);
            (*o) << std::endl;
            fr->print_contents(*o, export_alignment() ? '-' : 0x00);
            (*o) << std::endl;
        }
        (*o) << std::endl;
        if (!out) {
            delete o;
        }
    }
    return false;
}

const char* Output::name_impl() const {
    return "Output block set";
}

}

