/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>

#include "Fragment.hpp"
#include "Sequence.hpp"
#include "complement.hpp"

namespace bloomrepeats {

Fragment::Fragment(SequencePtr seq, size_t min_pos, size_t max_pos, int ori):
    seq_(seq), min_pos_(min_pos), max_pos_(max_pos), ori_(ori)
{ }

size_t Fragment::length() const {
    return max_pos() - min_pos() + 1;
}

size_t Fragment::begin_pos() const {
    return ori() == 1 ? min_pos() : max_pos();
}

const char* Fragment::begin() const {
    size_t l = length();
    return seq()->get(begin_pos(), l);
}

size_t Fragment::end_pos() const {
    return ori() == 1 ? max_pos() + 1 : min_pos() - 1;
}

const char* Fragment::end() const {
    return begin() + length() * ori();
}

std::ostream& operator<<(std::ostream& o, const Fragment& f) {
    for (const char* c = f.begin(); c != f.end(); c += f.ori()) {
        o << (f.ori() == 1 ? *c : complement(*c));
    }
    return o;
}

}

