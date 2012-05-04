/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SEQUENCE_HPP_
#define BR_SEQUENCE_HPP_

#include <string>

#include "global.hpp"

namespace bloomrepeats {

class Sequence {
public:
    virtual const char* get(size_t start, size_t& length) const = 0;
    virtual size_t approximate_size() const = 0;
};

class InMemorySequence : public Sequence {
public:
    // reads first sequence
    InMemorySequence(const std::string& filename);

    const char* get(size_t start, size_t& length) const;

    size_t approximate_size() const;

private:
    std::string data_;
};

}

#endif

