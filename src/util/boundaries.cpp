/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "boundaries.hpp"
#include "throw_assert.hpp"

namespace npge {

pos_t avg_element(const Boundaries& boundaries) {
    ASSERT_FALSE(boundaries.empty());
    pos_t result = 0;
    pos_t reminders = 0;
    BOOST_FOREACH (pos_t i, boundaries) {
        result += i / boundaries.size();
        reminders += i % boundaries.size();
    }
    result += reminders / boundaries.size();
    return result;
}

double avg_element_double(const Floats& floats) {
    ASSERT_FALSE(floats.empty());
    double sum = 0;
    BOOST_FOREACH (double f, floats) {
        sum += f;
    }
    return floats.size() ? sum / floats.size() : 0;
}

double avg_element_double(const Boundaries& boundaries) {
    ASSERT_FALSE(boundaries.empty());
    double sum = 0;
    BOOST_FOREACH (pos_t b, boundaries) {
        sum += b;
    }
    return boundaries.size() ? (sum / boundaries.size()) : 0;
}

Decimal avg_element_double(const Decimals& decimals) {
    Boundaries boundaries;
    BOOST_FOREACH (Decimal d, decimals) {
        boundaries.push_back(d.impl_);
    }
    int avg_impl = avg_element(boundaries);
    Decimal result;
    result.impl_ = avg_impl;
    return result;
}

pos_t nearest_element(const Boundaries& boundaries,
                      pos_t pos) {
    ASSERT_TRUE(boundaries.begin() != boundaries.end());
    Boundaries::const_iterator it = std::lower_bound(boundaries.begin(),
                                    boundaries.end(), pos);
    if (it == boundaries.end()) {
        // last
        --it;
        pos_t left = *it;
        return left;
    }
    pos_t right = *it;
    if (right == pos) {
        return pos;
    }
    if (it == boundaries.begin()) {
        return right;
    }
    --it;
    pos_t left = *it;
    if (right - pos <= pos - left) {
        return right;
    } else {
        return left;
    }
}

void select_boundaries(Boundaries& boundaries,
                       int min_distance,
                       pos_t length) {
    std::sort(boundaries.begin(), boundaries.end());
    Boundaries new_boundaries;
    Boundaries boundaries_nearby;
    pos_t prev = -1;
    BOOST_FOREACH (pos_t boundary, boundaries) {
        if (!boundaries_nearby.empty()) {
            if (boundary - boundaries_nearby[0] < min_distance) {
                boundaries_nearby.push_back(boundary);
            } else {
                new_boundaries.push_back(avg_element(boundaries_nearby));
                boundaries_nearby.clear();
                prev = boundary;
            }
        } else if (prev != -1 && boundary - prev < min_distance) {
            boundaries_nearby.push_back(prev);
            boundaries_nearby.push_back(boundary);
        } else {
            if (prev != -1) {
                new_boundaries.push_back(prev);
            }
            prev = boundary;
        }
    }
    if (!boundaries_nearby.empty()) {
        new_boundaries.push_back(avg_element(boundaries_nearby));
    } else if (prev != -1) {
        new_boundaries.push_back(prev);
    }
    if (new_boundaries.size() >= 1 && new_boundaries[0] - 0 < min_distance) {
        new_boundaries[0] = 0;
    }
    if (!new_boundaries.empty()) {
        ASSERT_LTE(new_boundaries.back(), length);
    }
    if (new_boundaries.size() >= 1 &&
            length - new_boundaries.back() < min_distance) {
        new_boundaries.back() = length;
    }
    boundaries.swap(new_boundaries);
}

}

