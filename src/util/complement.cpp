/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "complement.hpp"
#include "char_to_size.hpp"
#include "make_hash.hpp"

namespace npge {

void complement(std::string& str) {
    BOOST_FOREACH (char& c, str) {
        c = complement(c);
    }
    std::reverse(str.begin(), str.end());
}

hash_t complement_hash(hash_t hash, int letters_number) {
    hash_t result = 0;
    int digits = std::min(letters_number, MAX_ANCHOR_SIZE);
    // here we are using the fact that
    // complement_letter(x) = x ^ 1
    // r is number of full times our fragment fits in hash % 2
    size_t r = (letters_number / MAX_ANCHOR_SIZE) % 2;
    // position after which xor'ring rule changes
    size_t xor_swap_pos = letters_number % MAX_ANCHOR_SIZE;
    // xori'ring rule for positions smaller than xor_swap_pos
    size_t smal_xor = r ^ 1;
    // xori'ring rule for positions >= than xor_swap_pos
    size_t great_xor = r;
    for (int i = 0; i < digits; i++) {
        int dest = (letters_number - 1 - i) % MAX_ANCHOR_SIZE;
        hash_t last_2_bits = hash & 0x03;
        if (i < xor_swap_pos) {
            last_2_bits ^= smal_xor;
        } else {
            last_2_bits ^= great_xor;
        }
        result ^= last_2_bits << (2 * dest);
        // cut last 2 bits
        hash = hash >> 2;
    }
    return result;
}

}

