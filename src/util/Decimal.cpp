/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Decimal.hpp"

namespace npge {

std::ostream& operator<<(std::ostream& o, const Decimal& d) {
    o << d.to_s();
    return o;
}

}

