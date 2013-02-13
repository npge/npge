/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <vector>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "ImportBlastHits.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Exception.hpp"
#include "po.hpp"

namespace bloomrepeats {

ImportBlastHits::ImportBlastHits(const BlockSetPtr& block_set, int min_length,
                                 float min_ident, float max_evalue):
    min_length_(min_length), min_ident_(min_ident), max_evalue_(max_evalue) {
    set_other(block_set);
}

void ImportBlastHits::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("blast-hits", po::value<Files>()->required(),
     "results of blast -m 8")
    ("blast-min-length", po::value<int>()->default_value(min_length()),
     "min length of blast hit")
    ("blast-min-ident", po::value<float>()->default_value(min_ident()),
     "min ident of blast hit")
    ("blast-max-evalue", po::value<float>()->default_value(max_evalue()),
     "max e-value of blast hit")
   ;
}

void ImportBlastHits::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("blast-hits")) {
        set_input_files(vm["blast-hits"].as<Files>());
    }
    int min_length = vm["blast-min-length"].as<int>();
    if (min_length < 0) {
        throw Exception("'blast-min-length' must be >= 0");
    }
    set_min_length(min_length);
    float min_ident = vm["blast-min-ident"].as<float>();
    if (min_ident < 0 || min_ident > 1) {
        throw Exception("'blast-min-ident' must be in [0, 1]");
    }
    set_min_ident(min_ident);
    float max_evalue = vm["blast-max-evalue"].as<float>();
    if (max_evalue < 0) {
        throw Exception("'blast-max-evalue' must be >= 0");
    }
    set_max_evalue(max_evalue);
}

struct BlastItem {
    std::string id;
    int start;
    int stop;

    bool operator==(const BlastItem& o) const {
        return id == o.id && start == o.start && stop == o.stop;
    }

    bool operator!=(const BlastItem& o) const {
        return !(*this == o);
    }
};

struct BlastHit {
    BlastHit(std::string line) {
        using namespace boost::algorithm;
        trim(line);
        std::vector<std::string> parts;
        split(parts, line, is_any_of("\t"));
        items[0].id = parts[0];
        items[1].id = parts[1];
        ident = boost::lexical_cast<float>(parts[2]);
        length = boost::lexical_cast<int>(parts[3]);
        mismatches = boost::lexical_cast<int>(parts[4]);
        gap_openings = boost::lexical_cast<int>(parts[5]);
        items[0].start = boost::lexical_cast<int>(parts[6]);
        items[0].stop = boost::lexical_cast<int>(parts[7]);
        items[1].start = boost::lexical_cast<int>(parts[8]);
        items[1].stop = boost::lexical_cast<int>(parts[9]);
        evalue = boost::lexical_cast<double>(parts[10]);
        //bit_score = boost::lexical_cast<int>(parts[11]);
    }

    BlastItem items[2];
    float ident;
    int length;
    int mismatches;
    int gap_openings;
    double evalue;
    //int bit_score;
};

typedef std::map<std::string, Block*> NameToBlock;

static void add_blast_item(const BlockSet* bs, const NameToBlock& name2block,
                           Block* new_block, const BlastItem& item) {
    Fragment* f = bs->fragment_from_id(item.id);
    if (f) {
        new_block->insert(f->subfragment(item.start, item.stop));
        delete f;
    } else {
        NameToBlock::const_iterator it = name2block.find(item.id);
        BOOST_ASSERT(it != name2block.end());
        const Block* block = it->second;
        BOOST_ASSERT(block);
        BOOST_FOREACH (Fragment* fr, *block) {
            AlignmentRow* row = fr->row();
            BOOST_ASSERT(row);
            int start = row->nearest_in_fragment(item.start);
            BOOST_ASSERT(start != -1);
            int stop = row->nearest_in_fragment(item.stop);
            BOOST_ASSERT(stop != -1);
            new_block->insert(fr->subfragment(start, stop));
        }
    }
}

bool ImportBlastHits::run_impl() const {
    int size_before = block_set()->size();
    NameToBlock name2block;
    BOOST_FOREACH (Block* block, *other()) {
        name2block[block->name()] = block;
    }
    BlockSet* bs = other().get();
    BOOST_FOREACH (std::string file_name, input_files()) {
        std::ifstream input_file(file_name.c_str());
        for (std::string line; std::getline(input_file, line);) {
            BlastHit hit(line);
            if (hit.items[0] != hit.items[1] &&
                    hit.length >= min_length() &&
                    hit.ident >= min_ident() &&
                    hit.evalue <= max_evalue()) {
                Block* new_block = new Block;
                add_blast_item(bs, name2block, new_block, hit.items[0]);
                add_blast_item(bs, name2block, new_block, hit.items[1]);
                block_set()->insert(new_block);
            }
        }
    }
    return block_set()->size() > size_before;
}

const char* ImportBlastHits::name_impl() const {
    return "Import blast hits";
}

}

