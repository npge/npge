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

/** Convert char ('a', 't', 'g' or 'c') into size_t representation */
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

/** Make hash value from fragment of sequence.
\param start Beginning of the fragment
\param length Length of the fragment
\param ori Orientation of the fragment (1 or -1)
*/
inline size_t make_hash(const char* start, size_t length, int ori) {
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

#endif

