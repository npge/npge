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
#include <boost/algorithm/string/replace.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "AnchorFinder.hpp"
#include "JoinApprover.hpp"
#include "Exception.hpp"
#include "po.hpp"

using namespace bloomrepeats;

int main(int argc, char** argv) {
    po::options_description desc("Options");
    desc.add_options()
    ("help,h", "produce help message")
    ("input-file,i", po::value<std::vector<std::string> >()->required(),
     "input fasta file(s)")
    ("out-file,o", po::value<std::string>(), "output file with all blocks")
    ("out-mask", po::value<std::string>(),
     "mask of output files (${block} is replaced with block name)")
   ;
    po::positional_options_description pod;
    pod.add("input-file", -1);
    AnchorFinder anchor_finder;
    anchor_finder.add_options(desc);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(desc).positional(pod).run(), vm);
    } catch (std::exception& e) {
        std::cerr << argv[0] << ": error while parsing options: "
                  << std::endl << "  " << e.what() << std::endl;
        return 255;
    }
    if (vm.count("help")) {
        std::cout << "Usage:" << std::endl;
        std::cout << argv[0] << " [-i] input.fasta [options]" << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }
    try {
        po::notify(vm);
    } catch (std::exception& e) {
        std::cerr << argv[0] << ": error while notifying options: "
                  << std::endl << "  " << e.what() << std::endl;
        return 255;
    }
    try {
        anchor_finder.apply_options(vm);
    } catch (Exception& e) {
        std::cerr << argv[0] << ": error while setting up anchor finder: "
                  << std::endl << "  " << e.what() << std::endl;
        return 255;
    }
    BOOST_FOREACH (std::string file_name,
                  vm["input-file"].as<std::vector<std::string> >()) {
        std::ifstream input_file(file_name.c_str());
        while (true) {
            SequencePtr seq(new CompactSequence(input_file));
            if (seq->size() > 0) {
                anchor_finder.add_sequence(seq);
            } else {
                break;
            }
        }
    }
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    anchor_finder.set_block_set(block_set);
    anchor_finder.run();
    block_set->connect_fragments();
    block_set->resolve_overlaps();
    block_set->join(0);
    block_set->filter(10);
    block_set->expand_blocks_by_fragments();
    block_set->expand_blocks();
    block_set->filter(100);
    JoinApprover dist_1000(1000);
    block_set->join(&dist_1000);
#ifndef NDEBUG
    block_set->connect_fragments();
    BOOST_ASSERT(!block_set->overlaps());
#endif
    if (!vm["out-mask"].empty()) {
        std::string mask = vm["out-mask"].as<std::string>();
        BOOST_ASSERT(mask.find("${block}") != std::string::npos);
        BOOST_FOREACH (BlockPtr b, *block_set) {
            using namespace boost::algorithm;
            std::string path = replace_all_copy(mask, "${block}", b->name());
            std::ofstream o(path.c_str());
            o << *b << std::endl;
        }
    }
    if (!vm["out-file"].empty()) {
        std::string path = vm["out-file"].as<std::string>();
        std::ofstream o(path.c_str());
        o << *block_set << std::endl;
    }
    if (vm["out-file"].empty() && vm["out-mask"].empty()) {
        std::cout << *block_set << std::endl;
    }
}

