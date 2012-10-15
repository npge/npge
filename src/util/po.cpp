/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <exception>
#include <iostream>
#include <boost/foreach.hpp>

#include "po.hpp"

namespace bloomrepeats {

AddUniqueOptions::AddUniqueOptions(po::options_description& desc):
    po::options_description_easy_init(this),
    desc_(desc)
{ }

AddUniqueOptions::~AddUniqueOptions() {
    BOOST_FOREACH (boost::shared_ptr<po::option_description> opt, options()) {
        if (!desc_.find_nothrow(opt->long_name(), /* approx */ false)) {
            desc_.add(opt);
        }
    }
}

AddUniqueOptions add_unique_options(po::options_description& desc) {
    return AddUniqueOptions(desc);
}

void add_general_options(po::options_description& desc) {
    add_unique_options(desc)
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
    } catch (...) {
        std::cerr << argv[0] << ": error while parsing options: "
                  << std::endl << "  " << "Unknown error" << std::endl;
        return 255;
    }
    if (vm.count("help")) {
        std::cout << "Usage:" << std::endl;
        std::cout << argv[0] << " [options]";
        if (pod.max_total_count() > 0) {
            std::cout << " input";
        }
        std::cout << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }
    try {
        po::notify(vm);
    } catch (std::exception& e) {
        std::cerr << argv[0] << ": error while notifying options: "
                  << std::endl << "  " << e.what() << std::endl;
        return 255;
    } catch (...) {
        std::cerr << argv[0] << ": error while notifying options: "
                  << std::endl << "  " << "Unknown error" << std::endl;
        return 255;
    }
    return 0;
}

}

