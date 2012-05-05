/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_COMPLEMENT_HPP_
#define BR_COMPLEMENT_HPP_

namespace bloomrepeats {

inline char complement(char c) {
    if (c == 'a') {
        return 't';
    } else if (c == 't') {
        return 'a';
    } else if (c == 'g') {
        return 'c';
    } else { // if (c == 'c') {
        return 'g';
    }
}

}

#endif

