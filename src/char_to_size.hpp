/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CHAR_TO_SIZE_HPP_
#define BR_CHAR_TO_SIZE_HPP_

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

/** Convert size_t representation to char ('a', 't', 'g' or 'c') */
inline char size_to_char(size_t s) {
    if (s == 0) {
        return 'a';
    } else if (s == 1) {
        return 't';
    } else if (s == 2) {
        return 'g';
    } else { // if (c == 3) {
        return 'c';
    }
}

}

#endif

