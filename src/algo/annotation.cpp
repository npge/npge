/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/algorithm/string/predicate.hpp>

#include "annotation.hpp"

namespace npge {

bool is_id(const std::string& line) {
    using namespace boost::algorithm;
    return starts_with(line, "ID ") ||
        starts_with(line, "LOCUS ");
}

}
