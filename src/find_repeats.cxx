/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cassert>
#include <iostream>
#include <boost/lexical_cast.hpp>

#include "Sequence.hpp"
#include "BloomFilter.hpp"

using namespace bloomrepeats;

int main(int argc, char** argv) {
    assert(argc >= 2);
    size_t repeat_length = 18;
    if (argc >= 3) {
        repeat_length = boost::lexical_cast<int>(argv[2]);
    }
    InMemorySequence seq(argv[1]);
    float error_prob = 1.0 / seq.approximate_size();
    BloomFilter filter(seq.approximate_size(), error_prob);
    for (size_t start = 0; ; start++) {
        size_t length = repeat_length;
        const char* data = seq.get(start, length);
        if (length == repeat_length) {
            if (filter.test_and_add(data, length)) {
                std::cout << std::string(data, length) << std::endl;
            }
        } else {
            break;
        }
    }
}

