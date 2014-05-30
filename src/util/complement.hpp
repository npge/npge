/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_COMPLEMENT_HPP_
#define BR_COMPLEMENT_HPP_

#include <string>
#include <cctype>

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

}

#endif

