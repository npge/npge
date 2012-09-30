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
    if (c == 'a') {
        return 't';
    } else if (c == 't') {
        return 'a';
    } else if (c == 'g') {
        return 'c';
    } else if (c == 'c') {
        return 'g';
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

