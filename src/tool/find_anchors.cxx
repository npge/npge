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
    AnchorFinder anchor_finder;
    po::positional_options_description pod;
    pod.add("input-file", -1);
    BlockSet::add_output_options(desc);
    BlockSet::add_pangenome_options(desc);
    anchor_finder.add_options(desc);
    po::variables_map vm;
    int error = read_options(argc, argv, vm, desc, pod);
    if (error) {
        return error;
    }
    try {
        anchor_finder.apply_options(vm);
    } catch (Exception& e) {
        std::cerr << argv[0] << ": error while setting up anchor finder: "
                  << std::endl << "  " << e.what() << std::endl;
        return 255;
    }
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    std::vector<SequencePtr> seqs;
    Sequence::read_all_files(vm, seqs);
    block_set->add_sequences(seqs);
    anchor_finder.set_block_set(block_set);
    anchor_finder.run();
    block_set->make_output(vm);
}

