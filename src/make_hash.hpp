/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MAKE_HASH_HPP_
#define BR_MAKE_HASH_HPP_

#include "complement.hpp"

namespace bloomrepeats {

inline size_t char_to_size(char c) {
    if (c == 'a') {
        return 0;
    } else if (c == 't') {
        return 1;
    } else if (c == 'g') {
        return 2;
    } else { // if (c == 'c') {
        return 3;
    }
}

inline size_t make_hash(size_t hash_mul, const char* start,
                        size_t length, int ori) {
    size_t result = 1;
    const char* end = start + length * ori;
    for (const char* i = start; i != end; i += ori) {
        size_t value = char_to_size(ori == 1 ? *i : complement(*i));
        result *= hash_mul;
        result ^= value;
    }
    return result;
}

}

#endif

