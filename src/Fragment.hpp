/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FRAGMENT_HPP_
#define BR_FRAGMENT_HPP_

#include <iosfwd>
#include <string>
#include <boost/weak_ptr.hpp>

#include "global.hpp"

namespace bloomrepeats {

class Fragment {
public:
    Fragment();

    Fragment(SequencePtr seq, size_t min_pos, size_t max_pos, int ori);

    SequencePtr seq() const {
        return seq_;
    }

    BlockPtr block() const;

    size_t min_pos() const {
        return min_pos_;
    }

    void set_min_pos(size_t min_pos) {
        min_pos_ = min_pos;
    }

    size_t max_pos() const {
        return max_pos_;
    }

    void set_max_pos(size_t max_pos) {
        max_pos_ = max_pos;
    }

    int ori() const {
        return ori_;
    }

    size_t length() const;

    void set_ori(int ori) {
        ori_ = ori;
    }

    size_t begin_pos() const;

    const char* begin() const;

    size_t end_pos() const;

    const char* end() const;

    std::string str() const;

private:
    SequencePtr seq_;
    size_t min_pos_;
    size_t max_pos_;
    int ori_;
    boost::weak_ptr<Block> block_;
};

std::ostream& operator<<(std::ostream& o, const Fragment& fragment);

}

#endif

