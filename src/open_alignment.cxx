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

using namespace bloomrepeats;

int main(int argc, char** argv) {
    BOOST_ASSERT(argc >= 2);
    std::string filename = argv[1];
    std::ifstream input_file(filename.c_str());
    std::vector<SequencePtr> seqs;
    std::vector<FragmentPtr> frags;
    Alignment alignment;
    while (true) {
        SequencePtr seq(new InMemorySequence(input_file));
        if (seq->size() > 0) {
            // take first part
            seq->set_name(seq->name().substr(0, seq->name().find('_')));
            seqs.push_back(seq);
            FragmentPtr frag = new Fragment(seq, 0, seq->size());
            frags.push_back(frag);
            alignment.add_fragment(frag);
        } else {
            break;
        }
    }
    input_file.clear();
    input_file.seekg(0, std::ios::beg);
    input_file >> alignment;
    std::cout << alignment;
    BOOST_FOREACH (FragmentPtr frag, frags) {
        delete frag;
    }
    frags.clear();
}

