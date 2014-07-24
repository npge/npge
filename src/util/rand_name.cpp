/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <ctime>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "rand_name.hpp"

namespace npge {

int make_seed() {
    // http://stackoverflow.com/a/17623751
    using namespace boost::posix_time;
    ptime now = microsec_clock::local_time();
    int mcsec = now.time_of_day().total_milliseconds();
    int sec = time(NULL);
    return mcsec ^ sec;
}

static struct Srander {
    Srander() {
        std::srand(make_seed());
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

