/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "complement.hpp"
#include "char_to_size.hpp"

namespace npge {

void complement(std::string& str) {
    BOOST_FOREACH (char& c, str) {
        c = complement(c);
    }
    std::reverse(str.begin(), str.end());
}

hash_t complement_hash(hash_t hash, int letters_number) {
    hash_t result = 0;
    for (int i = 0; i < letters_number; i++) {
        int dest = letters_number - 1 - i;
        int last_2_bits = hash & 0x03;
        hash_t letter = complement_letter(last_2_bits);
        result |= letter << (2 * dest);
        // cut last 2 bits
        hash = hash >> 2;
    }
    return result;
}

}

