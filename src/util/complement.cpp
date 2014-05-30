/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "complement.hpp"

namespace npge {

void complement(std::string& str) {
    BOOST_FOREACH (char& c, str) {
        c = complement(c);
    }
    std::reverse(str.begin(), str.end());
}

}

