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
#include "Exception.hpp"
#include "po.hpp"

using namespace bloomrepeats;

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
    for (std::string line; std::getline(blast_hits_file, line);) {
        using namespace boost::algorithm;
        trim(line);
        std::vector<std::string> parts;
        split(parts, line, is_any_of("\t"));
        std::string f1_id = parts[0];
        std::string f2_id = parts[1];
        float ident = boost::lexical_cast<float>(parts[2]);
        int length = boost::lexical_cast<int>(parts[3]);
        int mismatches = boost::lexical_cast<int>(parts[4]);
        int gap_openings = boost::lexical_cast<int>(parts[5]);
        int f1_start = boost::lexical_cast<int>(parts[6]);
        int f1_end = boost::lexical_cast<int>(parts[7]);
        int f2_start = boost::lexical_cast<int>(parts[8]);
        int f2_end = boost::lexical_cast<int>(parts[9]);
        //double evalue = boost::lexical_cast<int>(parts[10]);
        //int bit_score = boost::lexical_cast<int>(parts[11]);
        FragmentPtr f1 = id2fragment[f1_id];
        FragmentPtr f2 = id2fragment[f2_id];
        if (f1 && f2 && *f1 != *f2) { // FIXME
            BOOST_ASSERT(f1);
            BOOST_ASSERT(f2);
            BlockPtr new_block = new Block;
            new_block->insert(f1->subfragment(f1_start, f1_end));
            new_block->insert(f2->subfragment(f2_start, f2_end));
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
    new_blocks->make_output(vm);
}

