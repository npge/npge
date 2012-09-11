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
#include "Fragment.hpp"
#include "BlockSet.hpp"

using namespace bloomrepeats;

int main(int argc, char** argv) {
    BOOST_ASSERT(argc >= 3);
    std::string sequences_filename = argv[1];
    std::string blocks_filename = argv[2];
    std::ifstream input_file(sequences_filename.c_str());
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
#ifndef NDEBUG
    size_t seq_size_summ = 0;
#endif
    while (true) {
        SequencePtr seq(new InMemorySequence(input_file));
#ifndef NDEBUG
        seq_size_summ += seq->size();
#endif
        if (seq->size() > 0) {
            block_set->add_sequence(seq);
        } else {
            break;
        }
    }
    std::ifstream blocks_file(blocks_filename.c_str());
    blocks_file >> *block_set;
    block_set->connect_fragments();
    BlockSetPtr all = block_set->rest();
#ifndef NDEBUG
    BOOST_ASSERT(!all->overlaps());
    all->connect_fragments();
    BOOST_ASSERT(!all->overlaps());
#endif
    BOOST_FOREACH (BlockPtr block, *block_set) {
        all->insert(block->clone());
    }
#ifndef NDEBUG
    BOOST_ASSERT(!all->overlaps());
    all->connect_fragments();
    BOOST_ASSERT(!all->overlaps());
    size_t fr_size_summ = 0;
    BOOST_FOREACH (BlockPtr block, *all) {
        BOOST_FOREACH (FragmentPtr f, *block) {
            fr_size_summ += f->length();
        }
    }
    BOOST_ASSERT(fr_size_summ == seq_size_summ);
#endif
    std::cout << *all << std::endl;
}

