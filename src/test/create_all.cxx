/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <vector>
#include <boost/foreach.hpp>

#include "Meta.hpp"
#include "Processor.hpp"
#include "global.hpp"

using namespace npge;

int main() {
    Meta meta;
    std::vector<SharedProcessor> processors;
    BOOST_FOREACH (std::string key, meta.keys()) {
        processors.push_back(meta.get(key));
    }
    std::cout << "Please check tmp files and press Enter" << "\n";
    std::cin.get();
}

