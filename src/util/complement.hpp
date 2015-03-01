/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_COMPLEMENT_HPP_
#define NPGE_COMPLEMENT_HPP_

#include <string>
#include <cctype>

#include "global.hpp"

namespace npge {

/** Return complementary nucleotide or input, if it is not in 'atgc' */
inline char complement(char c) {
    char C = toupper(c);
    if (C == 'A') {
        return 'T';
    } else if (C == 'T') {
        return 'A';
    } else if (C == 'G') {
        return 'C';
    } else if (C == 'C') {
        return 'G';
    } else {
        return c;
    }
}

/** Complement sequence in place.
Each of chars is complemented using complement(char),
and then the sequence is reversed.
*/
void complement(std::string& str);

/** Complement hash */
hash_t complement_hash(hash_t hash, int letters_number);

}

#endif

