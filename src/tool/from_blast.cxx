/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "AddSequences.hpp"
#include "AddBlocks.hpp"
#include "UniqueNames.hpp"
#include "Output.hpp"
#include "Alignment.hpp"
#include "Exception.hpp"
#include "po.hpp"

using namespace bloomrepeats;

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
        //evalue = boost::lexical_cast<int>(parts[10]);
        //bit_score = boost::lexical_cast<int>(parts[11]);
    }

    BlastItem items[2];
    float ident;
    int length;
    int mismatches;
    int gap_openings;
    //double evalue;
    //int bit_score;
};

typedef std::map<int, int> Int2Int;
typedef std::map<std::string, Int2Int> Frag2Map;
Frag2Map frag2map;
std::map<std::string, Fragment*> id2fragment;
std::map<std::string, Block*> name2block;

static void add_map(const BlastItem& item, const Alignment& alignment) {
    BOOST_FOREACH (Fragment* f, *alignment.block()) {
        int index = alignment.index_of(f);
        BOOST_ASSERT(index != -1);
        frag2map[f->id()][item.start] =
            alignment.nearest_in_fragment(index, item.start);
        frag2map[f->id()][item.stop] =
            alignment.nearest_in_fragment(index, item.stop);
    }
}

static void add_blast_item(Block* new_block, const BlastItem& item) {
    Fragment* f = id2fragment[item.id];
    if (f) {
        new_block->insert(f->subfragment(item.start, item.stop));
    } else {
        Block* block = name2block[item.id];
        BOOST_ASSERT(block);
        BOOST_FOREACH (Fragment* f, *block) {
            int start = frag2map[f->id()][item.start];
            int stop = frag2map[f->id()][item.stop];
            new_block->insert(f->subfragment(start, stop));
        }
    }
}

int main(int argc, char** argv) {
    po::options_description desc("Options");
    add_general_options(desc);
    AddSequences adder;
    adder.add_options(desc);
    AddBlocks blocks_adder(/* keep_alignment */ true);
    blocks_adder.add_options(desc);
    po::positional_options_description pod;
    Output output;
    output.add_options(desc);
    desc.add_options()
    ("blast-hits", po::value<std::string>()->required(),
     "input file with blast hits")
   ;
    po::variables_map vm;
    int error = read_options(argc, argv, vm, desc, pod);
    if (error) {
        return error;
    }
    try {
        adder.apply_options(vm);
        blocks_adder.apply_options(vm);
        output.apply_options(vm);
    } catch (Exception& e) {
        std::cerr << argv[0] << ": " << e.what() << std::endl;
        return 255;
    }
    BlockSetPtr pangenome = boost::make_shared<BlockSet>();
    adder.apply(pangenome);
    blocks_adder.apply(pangenome);
    std::ifstream blast_hits_file(vm["blast-hits"].as<std::string>().c_str());
    std::vector<BlastHit> blast_hits;
    for (std::string line; std::getline(blast_hits_file, line);) {
        BlastHit hit(line);
        if (hit.length >= 100 && hit.ident >= 0.95) {
            blast_hits.push_back(BlastHit(line));
        }
    }
    BOOST_FOREACH (Block* block, *pangenome) {
        BOOST_ASSERT(block->alignment());
        const Alignment& alignment = *block->alignment();
        const std::string& block_name = block->name();
        name2block[block_name] = block;
        BOOST_FOREACH (const BlastHit& hit, blast_hits) {
            if (hit.items[0].id == block_name) {
                add_map(hit.items[0], alignment);
            }
            if (hit.items[1].id == block_name) {
                add_map(hit.items[1], alignment);
            }
        }
    }
    BlockSetPtr new_blocks = boost::make_shared<BlockSet>();
    BOOST_FOREACH (Block* block, *pangenome) {
        BOOST_FOREACH (Fragment* f, *block) {
            id2fragment[f->id()] = f;
        }
    }
    BOOST_FOREACH (const BlastHit& hit, blast_hits) {
        if (hit.items[0] != hit.items[1]) {
            Block* new_block = new Block;
            add_blast_item(new_block, hit.items[0]);
            add_blast_item(new_block, hit.items[1]);
            new_blocks->insert(new_block);
        }
    }
    UniqueNames names;
    names.apply(new_blocks);
    output.apply(new_blocks);
}

