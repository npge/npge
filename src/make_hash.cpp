/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "make_hash.hpp"

namespace bloomrepeats {

size_t make_hash(const char* start, size_t length, int ori) {
    size_t result = 0;
    for (int j = 0; j < length; j++) {
        const char* i = start + j * ori;
        size_t value = char_to_size(ori == 1 ? *i : complement(*i));
        const int POS_BITS = 2;
        const int BYTE_BITS = 8;
        result ^= value << ((j * POS_BITS) % (sizeof(size_t) * BYTE_BITS));
    }
    return result;
}

}

