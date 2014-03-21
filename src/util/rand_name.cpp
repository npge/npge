/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <ctime>

#include "rand_name.hpp"

namespace bloomrepeats {

static struct Srander {
    Srander() {
        std::srand(time(NULL));
    }
} srander;

const char* const RAND_NAME_ABC = "0123456789abcdef";
const int RAND_NAME_ABC_SIZE = 16;

std::string rand_name(int size) {
    std::string result;
    result.resize(size);
    for (size_t i = 0; i < size; i++) {
        int r = rand() / (RAND_MAX / RAND_NAME_ABC_SIZE + 1);
        if (i == 0 && r < 10) {
            r = 10;
        }
        result[i] = RAND_NAME_ABC[r];
    }
    return result;
}

}

