/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>
#include <algorithm>
#include <boost/assert.hpp>

#include "Fragment.hpp"
#include "Sequence.hpp"
#include "complement.hpp"
#include "make_hash.hpp"

namespace bloomrepeats {

Fragment::Fragment(SequencePtr seq, size_t min_pos, size_t max_pos, int ori):
    seq_(seq), min_pos_(min_pos), max_pos_(max_pos), ori_(ori)
{ }

BlockPtr Fragment::block() const {
    return block_.lock();
}

FragmentPtr Fragment::prev() const {
    return prev_.lock();
}

FragmentPtr Fragment::next() const {
    return next_.lock();
}

FragmentPtr Fragment::neighbour(int ori) const {
    return ori == 1 ? next() : prev();
}

FragmentPtr Fragment::logical_neighbour(int ori) const {
    return neighbour(this->ori() * ori);
}

bool Fragment::is_neighbour(const Fragment& other) const {
    return prev().get() == &other || next().get() == &other;
}

FragmentPtr Fragment::another_neighbour(const Fragment& other) const {
    BOOST_ASSERT(is_neighbour(other));
    return prev().get() == &other ? next() : prev();
}

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

void Fragment::inverse() {
    set_ori(ori() == 1 ? -1 : 1);
}

std::string Fragment::str() const {
    std::string result;
    if (ori() == 1) {
        result.assign(begin(), length());
    } else {
        result.reserve(length());
        for (const char* c = begin(); c != end(); c--) {
            result += complement(*c);
        }
    }
    return result;
}

std::string Fragment::substr(int min, int max) const {
    min = min < 0 ? length() + min : min;
    max = max < 0 ? length() + max : max;
    int l = max - min + 1;
    std::string result;
    if (ori() == 1) {
        result.assign(begin() + min, l);
    } else {
        result.reserve(l);
        for (const char* c = begin() - min; c >= begin() - max; c--) {
            result += complement(*c);
        }
    }
    return result;
}

size_t Fragment::hash() const {
    return make_hash(begin(), length(), ori());
}

void Fragment::shift_end(int shift) {
    if (ori() == 1) {
        set_max_pos(max_pos() + shift);
    } else { /* if (ori() == -1) */
        set_min_pos(min_pos() - shift);
    }
}

int Fragment::max_shift_end() const {
    return ori() == 1 ? seq()->size() - max_pos() - 1 : min_pos();
}

bool Fragment::valid() const {
    return min_pos() <= max_pos() && max_pos() < seq()->size();
}

bool Fragment::operator==(const Fragment& other) const {
    return min_pos() == other.min_pos() && max_pos() == other.max_pos() &&
           ori() == other.ori() && seq() == other.seq();
}

bool Fragment::operator!=(const Fragment& other) const {
    return !(*this == other);
}

char Fragment::raw_at(int pos) const {
    char raw = *(begin() + ori() * pos);
    return ori() == 1 ? raw : complement(raw);
}

char Fragment::at(int pos) const {
    return raw_at(pos >= 0 ? pos : length() + pos);
}

void Fragment::connect(FragmentPtr first, FragmentPtr second) {
    if (first->next() != second) {
        if (first->next()) {
            first->next()->prev_.reset();
        }
        if (second->prev()) {
            second->prev()->next_.reset();
        }
    }
    first->next_ = second;
    second->prev_ = first;
}

bool Fragment::can_merge(FragmentPtr one, FragmentPtr another) {
    return one->seq() == another->seq() && one->ori() == another->ori() &&
           one->is_neighbour(*another);
}

FragmentPtr Fragment::merge(FragmentPtr one, FragmentPtr another) {
    BOOST_ASSERT(can_merge(one, another));
    if (another->next() == one) {
        std::swap(one, another);
    }
    FragmentPtr new_fragment = boost::make_shared<Fragment>(one->seq(),
                               std::min(one->min_pos(), another->min_pos()),
                               std::max(one->max_pos(), another->max_pos()),
                               one->ori());
    if (one->prev()) {
        connect(one->prev(), new_fragment);
    }
    if (another->next()) {
        connect(new_fragment, another->next());
    }
    return new_fragment;
}

void Fragment::disconnect() {
    if (next()) {
        next()->prev_.reset();
    }
    if (prev()) {
        prev()->next_.reset();
    }
    next_.reset();
    prev_.reset();
}

size_t Fragment::common_positions(const Fragment& other) {
    size_t result = 0;
    if (seq() == other.seq()) {
        size_t max_min = std::max(min_pos(), other.min_pos());
        size_t min_max = std::min(max_pos(), other.max_pos());
        if (max_min <= min_max) {
            result = min_max - max_min + 1;
        }
    }
    return result;
}

std::ostream& operator<<(std::ostream& o, const Fragment& f) {
    for (const char* c = f.begin(); c != f.end(); c += f.ori()) {
        o << (f.ori() == 1 ? *c : complement(*c));
    }
    return o;
}

}

