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
#include "AnchorFinder.hpp"
#include "Alignment.hpp"
#include "Exception.hpp"
#include "po.hpp"

using namespace bloomrepeats;

struct BlastHit {
    BlastHit(std::string line) {
        using namespace boost::algorithm;
        trim(line);
        std::vector<std::string> parts;
        split(parts, line, is_any_of("\t"));
        f1_id = parts[0];
        f2_id = parts[1];
        ident = boost::lexical_cast<float>(parts[2]);
        length = boost::lexical_cast<int>(parts[3]);
        mismatches = boost::lexical_cast<int>(parts[4]);
        gap_openings = boost::lexical_cast<int>(parts[5]);
        f1_start = boost::lexical_cast<int>(parts[6]);
        f1_end = boost::lexical_cast<int>(parts[7]);
        f2_start = boost::lexical_cast<int>(parts[8]);
        f2_end = boost::lexical_cast<int>(parts[9]);
        //evalue = boost::lexical_cast<int>(parts[10]);
        //bit_score = boost::lexical_cast<int>(parts[11]);
    }

    std::string f1_id;
    std::string f2_id;
    float ident;
    int length;
    int mismatches;
    int gap_openings;
    int f1_start;
    int f1_end;
    int f2_start;
    int f2_end;
    //double evalue;
    //int bit_score;
};

int main(int argc, char** argv) {
    po::options_description desc("Options");
    add_general_options(desc);
    Sequence::add_input_options(desc);
    po::positional_options_description pod;
    BlockSet::add_output_options(desc);
    desc.add_options()
    ("pangenome", po::value<std::string>()->required(),
     "input file with existing pangenome")
    ("blast-hits", po::value<std::string>()->required(),
     "input file with blast hits")
   ;
    po::variables_map vm;
    int error = read_options(argc, argv, vm, desc, pod);
    if (error) {
        return error;
    }
    BlockSetPtr pangenome = boost::make_shared<BlockSet>();
    std::vector<SequencePtr> seqs;
    Sequence::read_all_files(vm, seqs);
    pangenome->add_sequences(seqs);
    std::ifstream pangenome_file(vm["pangenome"].as<std::string>().c_str());
    pangenome_file >> *pangenome;
    std::map<std::string, FragmentPtr> id2fragment;
    BOOST_FOREACH (BlockPtr block, *pangenome) {
        BOOST_FOREACH (FragmentPtr f, *block) {
            id2fragment[f->id()] = f;
        }
    }
    std::ifstream blast_hits_file(vm["blast-hits"].as<std::string>().c_str());
    BlockSetPtr new_blocks = boost::make_shared<BlockSet>();
    std::vector<BlastHit> blast_hits;
    for (std::string line; std::getline(blast_hits_file, line);) {
        blast_hits.push_back(BlastHit(line));
    }
    BOOST_FOREACH (const BlastHit& hit, blast_hits) {
        FragmentPtr f1 = id2fragment[hit.f1_id];
        FragmentPtr f2 = id2fragment[hit.f2_id];
        if (f1 && f2 && *f1 != *f2) { // FIXME
            BOOST_ASSERT(f1);
            BOOST_ASSERT(f2);
            BlockPtr new_block = new Block;
            new_block->insert(f1->subfragment(hit.f1_start, hit.f1_end));
            new_block->insert(f2->subfragment(hit.f2_start, hit.f2_end));
            std::cout << *new_block;
            std::cout << std::endl;
            delete new_block;
        }
        //FragmentPtr f1_part = new Fragment(f1);
        //if (f1_start > f1_end) {
        //    f1_part->inverse();
        //    int tmp = f1_start;
        //    f1_start = f1_end;
        //    f1_end = tmp;
        //}
        //f1_part->set_begin_pos(f1_part->begin_pos() + f1_start);
        //FragmentPtr f2_part = new Fragment(f1);
    }
    new_blocks->set_unique_block_names();
    new_blocks->make_output(vm);
}

