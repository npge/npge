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
#include "Connector.hpp"
#include "Union.hpp"
#include "UniqueNames.hpp"
#include "OverlapsResolver.hpp"
#include "Rest.hpp"

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
    Connector connector;
    connector.apply(block_set);
    BlockSetPtr all = boost::make_shared<BlockSet>();
    Rest r(block_set);
    r.apply(all);
#ifndef NDEBUG
    OverlapsResolver resolver;
    resolver.set_block_set(all);
    BOOST_ASSERT(!resolver.overlaps());
    connector.apply(all);
    BOOST_ASSERT(!resolver.overlaps());
#endif
    Union u(block_set);
    u.apply(all);
#ifndef NDEBUG
    BOOST_ASSERT(!resolver.overlaps());
    connector.apply(all);
    BOOST_ASSERT(!resolver.overlaps());
    size_t fr_size_summ = 0;
    BOOST_FOREACH (Block* block, *all) {
        BOOST_FOREACH (Fragment* f, *block) {
            fr_size_summ += f->length();
        }
    }
    BOOST_ASSERT(fr_size_summ == seq_size_summ);
#endif
    UniqueNames names;
    names.apply(all);
    std::cout << *all << std::endl;
}

