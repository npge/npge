/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/utility/binary.hpp>

#include "make_hash.hpp"
#include "char_to_size.hpp"

namespace npge {

const int POS_BITS = 2;
const int BYTE_BITS = 8;
const size_t LAST_TWO_BITS = BOOST_BINARY(11);

static size_t shift_in_hash(int pos) {
    return ((pos * POS_BITS) % (sizeof(size_t) * BYTE_BITS));
}

size_t make_hash(const char* start, size_t length, int ori) {
    size_t result = 0;
    for (int j = 0; j < length; j++) {
        const char* i = start + j * ori;
        size_t value = char_to_size(ori == 1 ? *i : complement(*i));
        result ^= value << shift_in_hash(j);
    }
    return result;
}

size_t reuse_hash(size_t old_hash, size_t length,
                  char remove_char, char add_char, bool forward) {
    size_t remove = char_to_size(remove_char);
    old_hash ^= remove << shift_in_hash(forward ? 0 : length - 1);
    int occupied = std::min(POS_BITS * length, BYTE_BITS * sizeof(size_t));
    if (forward) {
        old_hash = (old_hash >> POS_BITS) |
                   ((old_hash & LAST_TWO_BITS) << (occupied - POS_BITS));
    } else {
        old_hash = (old_hash << POS_BITS) |
                   ((old_hash >> (occupied - POS_BITS)) & LAST_TWO_BITS);
    }
    size_t add = char_to_size(add_char);
    old_hash ^= add << shift_in_hash(forward ? length - 1 : 0);
    return old_hash;
}

}

