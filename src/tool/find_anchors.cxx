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
#include "AnchorFinder.hpp"
#include "Exception.hpp"
#include "po.hpp"

using namespace bloomrepeats;

void print_anchor(Block* block) {
    Fragment* fragment = *block->begin();
    std::cout << *fragment << std::endl;
    delete block;
}

int main(int argc, char** argv) {
    po::options_description desc("Options");
    desc.add_options()
    ("help,h", "produce help message")
    ("input-file,i", po::value<std::vector<std::string> >()->required(),
     "input fasta file(s)")
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
    anchor_finder.set_anchor_handler(print_anchor);
    anchor_finder.run();
}

