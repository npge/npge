/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_COMPLEMENT_HPP_
#define BR_COMPLEMENT_HPP_

#include <string>

namespace bloomrepeats {

/** Return complementary nucleotide or input, if it is not in 'atgc' */
inline char complement(char c) {
    if (c == 'A') {
        return 'T';
    } else if (c == 'T') {
        return 'A';
    } else if (c == 'G') {
        return 'C';
    } else if (c == 'C') {
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

