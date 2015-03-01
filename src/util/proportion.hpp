/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PROPORTION_HPP_
#define NPGE_PROPORTION_HPP_

#include "global.hpp"

namespace npge {

/** Return part1 / total1 * total2 */
inline pos_t proportion(pos_t part1, pos_t total1,
                        pos_t total2) {
    if (total1 == 0) {
        return 0;
    }
    double percentage = double(part1) / double(total1);
    return percentage * total2 + 0.00000001;
}

}

#endif

