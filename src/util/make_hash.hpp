/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_MAKE_HASH_HPP_
#define NPGE_MAKE_HASH_HPP_

#include <algorithm>
#include <boost/utility/binary.hpp>

#include "char_to_size.hpp"
#include "global.hpp"

namespace npge {

const int MAX_ANCHOR_SIZE = sizeof(hash_t) * 8 / 2;

const int POS_BITS = 2;
const int BYTE_BITS = 8;
const hash_t LAST_TWO_BITS = BOOST_BINARY(11);

inline size_t shift_in_hash(int pos) {
    return ((pos * POS_BITS) % (sizeof(hash_t) * BYTE_BITS));
}

template<typename F, int ori>
hash_t make_hash_base(F f, pos_t length) {
    hash_t result = 0;
    for (pos_t j = 0; j < length; j++) {
        char c = f(j);
        size_t s = char_to_size(c) & LAST_TWO_BITS;
        if (ori == -1 && c != 'N') {
            s = complement_letter(s);
        }
        hash_t value = s;
        result ^= value << shift_in_hash(j);
    }
    return result;
}

template<int ori>
struct FChar {
    const char* start_;

    FChar(const char* start):
        start_(start) {
    }

    char operator()(pos_t pos) {
        if (ori == 1) {
            return start_[pos];
        } else {
            int p = -int(pos);
            return start_[p];
        }
    }
};

/** Make hash value from fragment of sequence.
\param start Beginning of the fragment
\param length Length of the fragment
\param ori Orientation of the fragment (1 or -1)
*/
inline hash_t make_hash(const char* start, pos_t length,
                        int ori = 1) {
    if (ori == 1) {
        typedef FChar<1> F;
        F fchar((start));
        return make_hash_base<F, 1>(fchar, length);
    } else {
        typedef FChar < -1 > F;
        F fchar((start));
        return make_hash_base < F, -1 > (fchar, length);
    }
}

/** Make hash value from previous hash value (optimization).
\param old_hash Previous hash value
\param length Length of the fragment
\param forward If moved in according with the fragment's direction
\param remove_char Nucleotide, removed from the hash value.
\param add_char New nucleotide, added to the hash value.
Nucleotides add_char and remove_char should be pre-complement'ed, if needed.
*/
inline hash_t reuse_hash(hash_t old_hash, pos_t length,
                         char remove_char, char add_char,
                         bool forward = true) {
    hash_t remove = char_to_size(remove_char) & LAST_TWO_BITS;
    old_hash ^= remove << shift_in_hash(forward ?
                                        0 : length - 1);
    int occupied = std::min(int(POS_BITS * length),
                            int(BYTE_BITS * sizeof(hash_t)));
    if (forward) {
        old_hash = (old_hash >> POS_BITS) |
                   ((old_hash & LAST_TWO_BITS) <<
                    (occupied - POS_BITS));
    } else {
        old_hash = (old_hash << POS_BITS) |
                   ((old_hash >> (occupied - POS_BITS)) &
                    LAST_TWO_BITS);
    }
    hash_t add = char_to_size(add_char) & LAST_TWO_BITS;
    old_hash ^= add << shift_in_hash(forward ? length - 1 : 0);
    return old_hash;
}

}

#endif

