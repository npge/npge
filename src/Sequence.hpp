/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SEQUENCE_HPP_
#define BR_SEQUENCE_HPP_

#include <iosfwd>
#include <string>

#include "global.hpp"

namespace bloomrepeats {

class Sequence {
public:
    static const int FIRST_ORI = -1;

    virtual size_t approximate_size() const = 0;

    void make_first_fragment(Fragment& fragment, size_t fragment_size,
                             int only_ori = 1) const;

    bool next_fragment(Fragment& fragment) const;

    bool next_fragment_keeping_ori(Fragment& fragment) const;

protected:
    virtual const char* get(size_t start, size_t& length) const = 0;

private:
    friend class Fragment;
};

class InMemorySequence : public Sequence {
public:
    // reads first sequence
    InMemorySequence(const std::string& filename, int);

    InMemorySequence(std::istream& input);

    InMemorySequence(const std::string& data);

    size_t approximate_size() const;

protected:
    const char* get(size_t start, size_t& length) const;

private:
    std::string data_;

    void read_from_file(std::istream& input);
};

}

#endif

