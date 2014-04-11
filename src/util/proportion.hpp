/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PROPORTION_HPP_
#define BR_PROPORTION_HPP_

namespace bloomrepeats {

/** Return part1 / total1 * total2 */
inline int proportion(int part1, int total1, int total2) {
    if (total1 == 0) {
        return 0;
    }
    double percentage = double(part1) / double(total1);
    return percentage * total2 + 0.00000001;
}

}

#endif

