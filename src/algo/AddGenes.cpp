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
#include "throw_assert.hpp"
#include "global.hpp"

namespace bloomrepeats {

AddGenes::AddGenes():
    file_reader_(this, "in-genes", "input database files with genes") {
    add_opt("product", "Append product name to block name", true);
    declare_bs("target", "Blockset where genes are added");
}

void AddGenes::run_impl() const {
    BlockSet& bs = *block_set();
    std::map<std::string, Sequence*> ac2seq;
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        ac2seq[seq->ac()] = seq.get();
    }
    bool use_product = opt_value("product").as<bool>();
    BOOST_FOREACH (std::istream& input_file, file_reader_) {
        Sequence* seq = 0;
        Block* b = 0;
        Block* locus_tag_block = 0;
        for (std::string line; std::getline(input_file, line);) {
            using namespace boost::algorithm;
            if (starts_with(line, "AC")) {
                Strings parts;
                split(parts, line, isspace, token_compress_on);
                ASSERT_GTE(parts.size(), 2);
                std::string ac = parts[1];
                ASSERT_EQ(ac[ac.size() - 1], ';');
                ac.resize(ac.size() - 1);
                seq = ac2seq[ac];
                ASSERT_TRUE(seq);
                b = 0;
                locus_tag_block = 0;
            } else if (b && starts_with(line,
                                        "FT                   /locus_tag")) {
                Strings parts;
                split(parts, line, is_any_of("\""));
                const std::string& locus_tag = parts[1];
                b->set_name(locus_tag);
                locus_tag_block = b;
                b = 0;
            } else if (use_product && locus_tag_block &&
                       starts_with(line,
                                   "FT                   /product")) {
                Strings parts;
                split(parts, line, is_any_of("\""));
                const std::string& product = parts[1];
                std::string locus_tag = locus_tag_block->name();
                locus_tag += " " + product;
                locus_tag_block->set_name(locus_tag);
                locus_tag_block = 0;
            } else if (starts_with(line, "FT   CDS")) {
                ASSERT_TRUE(seq);
                Strings parts;
                split(parts, line, isspace, token_compress_on);
                ASSERT_GTE(parts.size(), 3);
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
                ASSERT_GT(coords.size(), 4);
                Strings boundaries;
                split(boundaries, coords, is_any_of(".,"),
                      token_compress_on);
                ASSERT_GTE(boundaries.size(), 2);
                ASSERT_EQ(boundaries.size() % 2, 0);
                b = new Block;
                bs.insert(b);
                for (int i = 0; i < boundaries.size() / 2; i++) {
                    std::string& min_pos_str = boundaries[i * 2];
                    std::string& max_pos_str = boundaries[i * 2 + 1];
                    // <1375315..1375356
                    if (!isdigit(min_pos_str[0])) {
                        min_pos_str = min_pos_str.substr(1);
                    }
                    if (!isdigit(max_pos_str[0])) {
                        max_pos_str = max_pos_str.substr(1);
                    }
                    int min_pos = boost::lexical_cast<int>(min_pos_str) - 1;
                    int max_pos = boost::lexical_cast<int>(max_pos_str) - 1;
                    Fragment* f = new Fragment(seq, min_pos, max_pos, ori);
                    b->insert(f);
                }
            }
        }
    }
}

const char* AddGenes::name_impl() const {
    return "Add genes from EBI genes description";
}

}

