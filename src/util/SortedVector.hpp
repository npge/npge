/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SORTED_VECTOR_HPP_
#define NPGE_SORTED_VECTOR_HPP_

#include <vector>
#include <ostream>
#include <algorithm>
#include <boost/foreach.hpp>

#include "throw_assert.hpp"

namespace npge {

/** Sorted array */
template <typename E>
class SortedVector : public std::vector<E> {
public:
    /** Typedef for vector */
    typedef typename std::vector<E> BaseStdVector;

    /** Constant iterator */
    typedef typename BaseStdVector::const_iterator const_iterator;

    /** Iterator */
    typedef typename BaseStdVector::iterator iterator;

    using BaseStdVector::begin;
    using BaseStdVector::end;
    using BaseStdVector::front;
    using BaseStdVector::back;
    using BaseStdVector::insert;
    using BaseStdVector::erase;
    using BaseStdVector::size;
    using BaseStdVector::capacity;
    using BaseStdVector::swap;
    using BaseStdVector::push_back;
    using BaseStdVector::pop_back;
    using BaseStdVector::at;

    /** Sort */
    void sort() {
        std::sort(begin(), end());
    }

    /** Remove duplicates (except first) from sorted vector */
    void unique() {
        erase(std::unique(begin(), end()), end());
        if (size() < capacity() / 2) {
            BaseStdVector copy((begin()), end());
            swap(copy);
        }
    }

    /** Remove all duplicates from sorted vector */
    void remove_multiple() {
        SortedVector<E> new_this;
        for (int i = 0; i < size(); i++) {
            bool prev = i - 1 >= 0 && at(i) == at(i - 1);
            bool next = i + 1 < size() && at(i) == at(i + 1);
            if (!prev && !next) {
                new_this.push_back(at(i));
            }
        }
        SortedVector<E>::swap(new_this);
        ASSERT_TRUE(is_sorted_unique());
    }

    /** Sort and remove duplicates */
    void sort_unique() {
        sort();
        unique();
    }

    /** Return if it is sorted and uniqued */
    bool is_sorted_unique() const {
        if (size() <= 1) {
            return true;
        }
        const_iterator i1, i2;
        i1 = begin();
        i2 = i1;
        ++i2;
        for (; i2 < end(); ++i1, ++i2) {
            if (!(*i1 < *i2)) {
                return false;
            }
        }
        return true;
    }

    /** Shortcut for std::lower_bound */
    const_iterator lower_bound(const E& elem) const {
        return std::lower_bound(begin(), end(), elem);
    }

    /** Shortcut for std::upper_bound */
    const_iterator upper_bound(const E& elem) const {
        return std::upper_bound(begin(), end(), elem);
    }

    /** Return if the vector contains the element */
    bool has_elem(const E& e) const {
        const_iterator it = lower_bound(e);
        return it != end() && *it == e;
    }

    /** Add all elements from other vector to end */
    void extend(const SortedVector& other) {
        insert(end(), other.begin(), other.end());
    }
};

/** Streaming operator */
template<typename E>
std::ostream& operator<<(std::ostream& o, const SortedVector<E>& v) {
    BOOST_FOREACH (const E& elem, v) {
        o << elem << std::endl;
    }
    return o;
}

}

#endif

