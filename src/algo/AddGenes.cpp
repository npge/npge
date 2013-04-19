/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <istream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "AddGenes.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

void AddGenes::add_options_impl(po::options_description& desc) const {
    bloomrepeats::add_unique_options(desc)
    ("in-genes", po::value<Files>()->multitoken()->required(),
     "input database files with genes")
   ;
}

void AddGenes::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("in-genes")) {
        set_input_files(vm["in-genes"].as<Files>());
    }
}

bool AddGenes::run_impl() const {
    BlockSet& bs = *block_set();
    int size_before = bs.size();
    BOOST_FOREACH (std::istream& input_file, *this) {
        SequencePtr seq;
        for (std::string line; std::getline(input_file, line);) {
            using namespace boost::algorithm;
            if (starts_with(line, "AC")) {
                std::vector<std::string> parts;
                split(parts, line, isspace, token_compress_on);
                BOOST_ASSERT(parts.size() >= 2);
                std::string ac = parts[1];
                BOOST_ASSERT(ac[ac.size() - 1] == ';');
                ac.resize(ac.size() - 1);
                seq = bs.seq_from_name(ac);
            } else if (starts_with(line, "FT   gene")) {
                BOOST_ASSERT(seq);
                std::vector<std::string> parts;
                split(parts, line, isspace, token_compress_on);
                BOOST_ASSERT(parts.size() >= 3);
                std::string& coords = parts[2];
                int ori = 1;
                if (starts_with(coords, "complement(")) {
                    ori = -1;
                    int slice_begin = 11;
                    int slice_end = coords.size() - 1;
                    int slice_length = slice_end - slice_begin;
                    coords = coords.substr(slice_begin, slice_length);
                }
                BOOST_ASSERT(coords.size() > 4);
                if (coords[0] == '<') {
                    coords = coords.substr(1); // TODO
                }
                std::vector<std::string> min_max;
                split(min_max, coords, is_any_of("."), token_compress_on);
                BOOST_ASSERT(min_max.size() == 2);
                int min_pos = boost::lexical_cast<int>(min_max[0]);
                int max_pos = boost::lexical_cast<int>(min_max[1]);
                Fragment* f = new Fragment(seq, min_pos, max_pos, ori);
                Block* b = new Block;
                b->insert(f);
                bs.insert(b);
            }
        }
    }
    return bs.size() > size_before;
}

const char* AddGenes::name_impl() const {
    return "Add genes from EBI genes description";
}

}

