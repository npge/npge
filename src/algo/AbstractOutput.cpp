/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "AbstractOutput.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

static bool file_and_mask_check(AbstractOutput* p,
        std::string& message) {
    if (p->opt_value("file").as<std::string>() != "" &&
            p->opt_value("mask").as<std::string>() != "") {
        message = "both '" + p->opt_prefixed("file") +
                  "' and '" + p->opt_prefixed("mask") +
                  "' were specified";
        return false;
    } else {
        return true;
    }
}

static bool mask_check(AbstractOutput* p, std::string& message) {
    std::string mask = p->opt_value("mask").as<std::string>();
    if (mask != "" && mask.find("${block}") == std::string::npos) {
        message = "'" + p->opt_prefixed("mask") +
            "' must contain '${block}'";
        return false;
    } else {
        return true;
    }
}

AbstractOutput::AbstractOutput() {
    add_opt("file", "output file with all blocks", std::string());
    add_opt("mask", "mask of output files (${block} is "
            "replaced with block name)", std::string());
    add_opt_check(boost::bind(file_and_mask_check, this, _1));
    add_opt_check(boost::bind(mask_check, this, _1));
}

static struct BlockCompareName2 {
    bool operator()(const Block* b1, const Block* b2) const {
        typedef boost::tuple<int, const std::string&> Tie;
        return Tie(-b1->size(), b1->name()) < Tie(-b2->size(), b2->name());
    }
} bcn2;

bool AbstractOutput::run_impl() const {
    std::string file = opt_value("file").as<std::string>();
    std::string mask = opt_value("mask").as<std::string>();
    prepare();
    std::vector<Block*> blocks(block_set()->begin(), block_set()->end());
    if (mask.empty()) {
        std::sort(blocks.begin(), blocks.end(), bcn2);
    }
    boost::shared_ptr<std::ostream> out;
    if (mask.empty()) {
        out = name_to_ostream(file);
    }
    if (out) {
        print_header(*out);
    }
    BOOST_FOREACH (Block* b, blocks) {
        boost::shared_ptr<std::ostream> o = out;
        std::string path;
        if (!out) {
            using namespace boost::algorithm;
            path = replace_all_copy(mask, "${block}", b->name());
            o = name_to_ostream(path);
            print_header(*o);
        }
        BOOST_ASSERT(o);
        print_block(*o, b);
        if (!out) {
            print_footer(*o);
        }
    }
    if (out) {
        print_footer(*out);
    }
    return false;
}

void AbstractOutput::prepare() const
{ }

void AbstractOutput::print_header(std::ostream& o) const
{ }

void AbstractOutput::print_footer(std::ostream& o) const
{ }

}

