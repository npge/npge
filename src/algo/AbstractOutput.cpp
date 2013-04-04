/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "AbstractOutput.hpp"
#include "Exception.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

void AbstractOutput::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("file", po::value<std::string>()->default_value(file()),
     "output file with all blocks")
    ("mask", po::value<std::string>()->default_value(mask()),
     "mask of output files (${block} is replaced with block name)")
   ;
}

void AbstractOutput::apply_options_impl(const po::variables_map& vm) {
    if (vm.count(prefixed("file"))) {
        set_file(vm[prefixed("file")].as<std::string>());
    }
    if (vm.count(prefixed("mask"))) {
        set_mask(vm[prefixed("mask")].as<std::string>());
    }
    if (!file().empty() && !mask().empty()) {
        throw Exception("both '" + prefixed("file") +
                        "' and '" + prefixed("mask") + "' were specified");
    }
    if (!mask().empty() && mask().find("${block}") == std::string::npos) {
        throw Exception("'" + prefixed("mask") + "' must contain '${block}'");
    }
}

static struct BlockCompareName2 {
    bool operator()(const Block* b1, const Block* b2) const {
        typedef boost::tuple<int, const std::string&> Tie;
        return Tie(-b1->size(), b1->name()) < Tie(-b2->size(), b2->name());
    }
} bcn2;

bool AbstractOutput::run_impl() const {
    std::vector<Block*> blocks(block_set()->begin(), block_set()->end());
    if (mask().empty()) {
        std::sort(blocks.begin(), blocks.end(), bcn2);
    }
    std::ostream* out = mask().empty() ? name_to_ostream(file()).get() : 0;
    if (out) {
        print_header(*out);
    }
    BOOST_FOREACH (Block* b, blocks) {
        std::ostream* o = out;
        std::string path;
        if (!out) {
            using namespace boost::algorithm;
            path = replace_all_copy(mask(), "${block}", b->name());
            o = name_to_ostream(path).get();
            print_header(*o);
        }
        BOOST_ASSERT(o);
        print_block(*o, b);
        if (!out) {
            print_footer(*o);
            remove_ostream(path, /* remove_file */ false);
        }
    }
    if (out) {
        print_footer(*out);
        remove_ostream(file(), /* remove_file */ false);
    }
    return false;
}

void AbstractOutput::print_header(std::ostream& o) const
{ }

void AbstractOutput::print_footer(std::ostream& o) const
{ }

}

