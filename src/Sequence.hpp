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

    size_t size() const {
        return size_;
    }

    void make_first_fragment(Fragment& fragment, size_t fragment_size,
                             int only_ori = 1) const;

    bool next_fragment(Fragment& fragment) const;

    bool next_fragment_keeping_ori(Fragment& fragment) const;

    static void to_atgc(std::string& data);

protected:
    virtual char char_at(size_t index) const = 0;

    void set_size(size_t size) {
        size_ = size;
    }

private:
    size_t size_;

    friend class Fragment;
};

class InMemorySequence : public Sequence {
public:
    // reads first sequence
    InMemorySequence(const std::string& filename, int);

    InMemorySequence(std::istream& input);

    InMemorySequence(const std::string& data);

protected:
    virtual char char_at(size_t index) const;

private:
    std::string data_;

    void read_from_file(std::istream& input);
};

}

#endif

