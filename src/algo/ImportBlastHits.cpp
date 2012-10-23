/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

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
#include "Alignment.hpp"
#include "Exception.hpp"
#include "po.hpp"

namespace bloomrepeats {

ImportBlastHits::ImportBlastHits(const BlockSetPtr& block_set,
                                 int min_length, float min_ident):
    OtherBlockSet(block_set),
    min_length_(min_length), min_ident_(min_ident)
{ }

void ImportBlastHits::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("blast-hits", po::value<Files>()->required(),
     "results of blast -m 8")
    ("blast-min-length", po::value<int>()->default_value(min_length()),
     "min length of blast hit")
    ("blast-min-ident", po::value<float>()->default_value(min_ident()),
     "min ident of blast hit")
   ;
}

void ImportBlastHits::apply_options_impl(const po::variables_map& vm) {
    set_files(vm["blast-hits"].as<Files>());
    int min_length = vm["blast-min-length"].as<int>();
    if (min_length < 0) {
        throw Exception("blast-min-length' must be >= 0");
    }
    set_min_length(min_length);
    float min_ident = vm["blast-min-ident"].as<float>();
    if (min_ident < 0 || min_ident > 1) {
        throw Exception("blast-min-ident' must be in [0, 1]");
    }
    set_min_ident(min_ident);
}

struct ImportBlastHits::BlastItem {
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
        //evalue = boost::lexical_cast<int>(parts[10]);
        //bit_score = boost::lexical_cast<int>(parts[11]);
    }

    ImportBlastHits::BlastItem items[2];
    float ident;
    int length;
    int mismatches;
    int gap_openings;
    //double evalue;
    //int bit_score;
};

bool ImportBlastHits::run_impl() const {
    int size_before = block_set()->size();
    std::vector<BlastHit> blast_hits;
    BOOST_FOREACH (std::string file_name, files()) {
        std::ifstream input_file(file_name.c_str());
        for (std::string line; std::getline(input_file, line);) {
            BlastHit hit(line);
            if (hit.length >= min_length() && hit.ident >= min_ident()) {
                blast_hits.push_back(BlastHit(line));
            }
        }
    }
    BOOST_FOREACH (Block* block, *other()) {
        BOOST_ASSERT(block->alignment());
        const Alignment& alignment = *block->alignment();
        const std::string& block_name = block->name();
        name2block_[block_name] = block;
        BOOST_FOREACH (const BlastHit& hit, blast_hits) {
            if (hit.items[0].id == block_name) {
                add_map(hit.items[0], alignment);
            }
            if (hit.items[1].id == block_name) {
                add_map(hit.items[1], alignment);
            }
        }
    }
    BOOST_FOREACH (Block* block, *other()) {
        BOOST_FOREACH (Fragment* f, *block) {
            id2fragment_[f->id()] = f;
        }
    }
    BOOST_FOREACH (const BlastHit& hit, blast_hits) {
        if (hit.items[0] != hit.items[1]) {
            Block* new_block = new Block;
            add_blast_item(new_block, hit.items[0]);
            add_blast_item(new_block, hit.items[1]);
            block_set()->insert(new_block);
        }
    }
    frag2map_.clear();
    id2fragment_.clear();
    name2block_.clear();
    return block_set()->size() > size_before;
}

const char* ImportBlastHits::name_impl() const {
    return "Import blast hits";
}

void ImportBlastHits::add_map(const ImportBlastHits::BlastItem& item,
                              const Alignment& alignment) const {
    BOOST_FOREACH (Fragment* f, *alignment.block()) {
        int index = alignment.index_of(f);
        BOOST_ASSERT(index != -1);
        frag2map_[f->id()][item.start] =
            alignment.nearest_in_fragment(index, item.start);
        frag2map_[f->id()][item.stop] =
            alignment.nearest_in_fragment(index, item.stop);
    }
}

void ImportBlastHits::add_blast_item(Block* new_block,
                                     const ImportBlastHits::BlastItem&
                                     item) const {
    Fragment* f = id2fragment_[item.id];
    if (f) {
        new_block->insert(f->subfragment(item.start, item.stop));
    } else {
        Block* block = name2block_[item.id];
        BOOST_ASSERT(block);
        BOOST_FOREACH (Fragment* f, *block) {
            int start = frag2map_[f->id()][item.start];
            int stop = frag2map_[f->id()][item.stop];
            new_block->insert(f->subfragment(start, stop));
        }
    }
}

}

