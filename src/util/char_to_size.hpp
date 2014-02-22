/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CHAR_TO_SIZE_HPP_
#define BR_CHAR_TO_SIZE_HPP_

namespace bloomrepeats {

/** Number letters */
const int LETTERS_NUMBER = 5;

/** Letters */
enum Letter {
    A, T, C, G, N
};

/** Convert char ('A', 'T', 'G', 'C', 'N') into size_t representation.
Max value returned is LETTERS_NUMBER - 1.
*/
inline size_t char_to_size(char c) {
    if (c == 'A') {
        return 0;
    } else if (c == 'T') {
        return 1;
    } else if (c == 'G') {
        return 2;
    } else if (c == 'C') {
        return 3;
    } else { // if (c == 'N') {
        return 4;
    }
}

/** Convert size_t representation to char ('A', 'T', 'G', 'C', 'N') */
inline char size_to_char(size_t s) {
    if (s == 0) {
        return 'A';
    } else if (s == 1) {
        return 'T';
    } else if (s == 2) {
        return 'G';
    } else if (s == 3) {
        return 'C';
    } else { // if (s == 4) {
        return 'N';
    }
}

}

#endif

