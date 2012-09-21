/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <exception>
#include <iostream>

#include "po.hpp"

namespace bloomrepeats {

void add_general_options(po::options_description& desc) {
    desc.add_options()
    ("help,h", "produce help message")
   ;
}

int read_options(int argc, char** argv, po::variables_map& vm,
                 const po::options_description& desc,
                 const po::positional_options_description& pod) {
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
    return 0;
}

}

