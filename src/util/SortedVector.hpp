/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SORTED_VECTOR_HPP_
#define BR_SORTED_VECTOR_HPP_

#include <vector>
#include <algorithm>

namespace bloomrepeats {

/** A graph represented as sorted array of pairs */
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
    using BaseStdVector::erase;
    using BaseStdVector::size;
    using BaseStdVector::push_back;
    using BaseStdVector::pop_back;

    /** Sort and remove duplicates */
    void sort_unique() {
        std::sort(begin(), end());
        erase(std::unique(begin(), end()), end());
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
            if (*i1 >= *i2) {
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
};

}

#endif

