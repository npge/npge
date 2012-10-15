/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <exception>

#include "process.hpp"
#include "BlockSet.hpp"
#include "po.hpp"

namespace bloomrepeats {

int process(int argc, char** argv,
            ProcessorPtr processor,
            const std::string& name,
            const std::string& positional) {
    po::options_description desc(name);
    add_general_options(desc);
    processor->add_options(desc);
    po::positional_options_description pod;
    if (!positional.empty()) {
        pod.add(positional.c_str(), -1);
    }
    po::variables_map vm;
    int error = read_options(argc, argv, vm, desc, pod);
    if (error) {
        return error;
    }
    try {
        processor->apply_options(vm);
    } catch (std::exception& e) {
        std::cerr << argv[0] << ": error while applying options" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return 255;
    } catch (...) {
        std::cerr << argv[0] << ": error while applying options" << std::endl;
        std::cerr << "  Unknown error" << std::endl;
        return 255;
    }
    try {
        BlockSetPtr block_set = boost::make_shared<BlockSet>();
        processor->apply(block_set);
    } catch (std::exception& e) {
        std::cerr << argv[0] << ": algorithm error" << std::endl;
        std::cerr << "  " << e.what() << std::endl;
        return 255;
    } catch (...) {
        std::cerr << argv[0] << ": error while applying options" << std::endl;
        std::cerr << "  Unknown error" << std::endl;
        return 255;
    }
    return 0;
}

}

