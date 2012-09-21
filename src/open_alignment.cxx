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

#include "Alignment.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

using namespace bloomrepeats;

int main(int argc, char** argv) {
    BOOST_ASSERT(argc >= 3);
    std::string sequences_filename = argv[1];
    std::string blocks_filename = argv[2];
    std::ifstream input_file(sequences_filename.c_str());
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    while (true) {
        SequencePtr seq(new InMemorySequence(input_file));
        if (seq->size() > 0) {
            block_set->add_sequence(seq);
        } else {
            break;
        }
    }
    std::ifstream blocks_file(blocks_filename.c_str());
    blocks_file >> *block_set;
    Alignment alignment;
    BOOST_ASSERT(block_set->size() == 1);
    BOOST_FOREACH (FragmentPtr f, *block_set->front()) {
        alignment.add_fragment(f);
    }
    blocks_file.clear();
    blocks_file.seekg(0, std::ios::beg);
    blocks_file >> alignment;
    std::cout << alignment;
}

