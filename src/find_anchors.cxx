/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cassert>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>

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
    SequencePtr seq = boost::make_shared<InMemorySequence>(argv[1]);
    AnchorFinder anchor_finder;
    anchor_finder.add_sequnce(seq);
    anchor_finder.set_anchor_handler(print_anchor);
    anchor_finder.set_anchor_size(repeat_length);
    anchor_finder.run();
}

