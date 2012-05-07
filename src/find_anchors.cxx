/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cassert>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "AnchorFinder.hpp"

using namespace bloomrepeats;

void print_anchor(BlockPtr block) {
    FragmentPtr fragment = *block->begin();
    std::cout << *fragment << std::endl;
}

int main(int argc, char** argv) {
    assert(argc >= 2);
    size_t repeat_length = 20;
    if (argc >= 3) {
        repeat_length = boost::lexical_cast<int>(argv[2]);
    }
    int workers = 1;
    if (argc >= 4) {
        workers = boost::lexical_cast<int>(argv[3]);
    }
    AnchorFinder anchor_finder;
    std::ifstream input_file(argv[1]);
    while (true) {
        SequencePtr seq(new InMemorySequence(input_file));
        if (seq->size() > 0) {
            anchor_finder.add_sequence(seq);
        } else {
            break;
        }
    }
    anchor_finder.set_anchor_handler(print_anchor);
    anchor_finder.set_anchor_size(repeat_length);
    anchor_finder.set_workers(workers);
    anchor_finder.run();
}

