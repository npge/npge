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
    block_set->connect_fragments();
    BlockSetPtr rest = block_set->rest();
#ifndef NDEBUG
    BOOST_ASSERT(!rest->overlaps());
    rest->connect_fragments();
    BOOST_ASSERT(!rest->overlaps());
#endif
    BOOST_FOREACH (BlockPtr block, *block_set) {
        rest->insert(block->clone());
    }
#ifndef NDEBUG
    BOOST_ASSERT(!rest->overlaps());
    rest->connect_fragments();
    BOOST_ASSERT(!rest->overlaps());
#endif
    std::cout << *rest << std::endl;
}

