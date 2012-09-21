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
    desc.add_options()
    ("help,h", "produce help message")
   ;
    Sequence::add_input_options(desc);
    po::positional_options_description pod;
    pod.add("input-file", -1);
    AnchorFinder anchor_finder;
    anchor_finder.add_options(desc);
    BlockSet::add_output_options(desc);
    BlockSet::add_pangenome_options(desc);
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
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    anchor_finder.set_block_set(block_set);
    std::vector<SequencePtr> seqs;
    Sequence::read_all_files(vm, seqs);
    anchor_finder.add_sequences(seqs);
    anchor_finder.run();
    block_set->make_pangenome(vm);
#ifndef NDEBUG
    block_set->connect_fragments();
    BOOST_ASSERT(!block_set->overlaps());
#endif
    block_set->make_output(vm);
}

