/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <istream>
#include <map>
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
    std::map<std::string, Sequence*> ac2seq;
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        ac2seq[seq->ac()] = seq.get();
    }
    int size_before = bs.size();
    BOOST_FOREACH (std::istream& input_file, *this) {
        Sequence* seq = 0;
        Block* b = 0;
        for (std::string line; std::getline(input_file, line);) {
            using namespace boost::algorithm;
            if (starts_with(line, "AC")) {
                std::vector<std::string> parts;
                split(parts, line, isspace, token_compress_on);
                BOOST_ASSERT(parts.size() >= 2);
                std::string ac = parts[1];
                BOOST_ASSERT(ac[ac.size() - 1] == ';');
                ac.resize(ac.size() - 1);
                seq = ac2seq[ac];
                BOOST_ASSERT(seq);
                b = 0;
            } else if (starts_with(line, "FT                   /locus_tag")) {
                std::vector<std::string> parts;
                split(parts, line, is_any_of("\""));
                const std::string& locus_tag = parts[1];
                if (b) {
                    b->set_name(locus_tag);
                    b = 0;
                }
            } else if (starts_with(line, "FT   CDS")) {
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
                if (starts_with(coords, "join(")) {
                    int slice_begin = 5;
                    int slice_end = coords.size() - 1;
                    int slice_length = slice_end - slice_begin;
                    coords = coords.substr(slice_begin, slice_length);
                }
                BOOST_ASSERT(coords.size() > 4);
                std::vector<std::string> min_max;
                split(min_max, coords, is_any_of("."), token_compress_on);
                BOOST_ASSERT(min_max.size() >= 2);
                int last = min_max.size() - 1; // join(1..20,23..40)
                std::string& min_pos_str = min_max[0];
                std::string& max_pos_str = min_max[last];
                if (!isdigit(min_pos_str[0])) {
                    min_pos_str = min_pos_str.substr(1);
                }
                if (!isdigit(max_pos_str[0])) {
                    max_pos_str = max_pos_str.substr(1);
                }
                int min_pos = boost::lexical_cast<int>(min_pos_str) - 1;
                int max_pos = boost::lexical_cast<int>(max_pos_str) - 1;
                Fragment* f = new Fragment(seq, min_pos, max_pos, ori);
                b = new Block;
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

