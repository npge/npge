/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "StickBoundaries.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

StickBoundaries::StickBoundaries(int min_distance):
    min_distance_(min_distance)
{ }

void StickBoundaries::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("min-distance", po::value<int>()->default_value(min_distance()),
     "Min distance between fragment boundaries")
   ;
}

void StickBoundaries::apply_options_impl(const po::variables_map& vm) {
    int min_distance = vm["min-distance"].as<int>();
    if (min_distance < 0) {
        throw Exception("'min-distance' must be >= 0");
    }
    set_min_distance(min_distance);
}

typedef std::vector<size_t> Boundaries;

static size_t avg_element(const Boundaries& boundaries) {
    size_t result = 0;
    size_t reminders = 0;
    BOOST_FOREACH (size_t i, boundaries) {
        result += i / boundaries.size();
        reminders += i % boundaries.size();
    }
    result += reminders / boundaries.size();
    return result;
}

static size_t nearest_element(const Boundaries& boundaries, size_t pos) {
    BOOST_ASSERT(boundaries.begin() != boundaries.end());
    Boundaries::const_iterator it = std::lower_bound(boundaries.begin(),
                                    boundaries.end(), pos);
    if (it == boundaries.end()) {
        // last
        --it;
        size_t left = *it;
        return left;
    }
    size_t right = *it;
    if (right == pos) {
        return pos;
    }
    if (it == boundaries.begin()) {
        return right;
    }
    --it;
    size_t left = *it;
    if (right - pos < pos - left) {
        return right;
    } else {
        return left;
    }
}

static void select_boundaries(Boundaries& boundaries, int min_distance) {
    std::sort(boundaries.begin(), boundaries.end());
    Boundaries new_boundaries;
    Boundaries boundaries_nearby;
    size_t prev = -1;
    BOOST_FOREACH (size_t boundary, boundaries) {
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
    boundaries.swap(new_boundaries);
}

bool StickBoundaries::run_impl() const {
    typedef std::map<const Sequence*, Boundaries> Seq2Boundaries;
    Seq2Boundaries seq2b;
    BOOST_FOREACH (Block* block, *block_set()) {
        BOOST_FOREACH (Fragment* fragment, *block) {
            seq2b[fragment->seq()].push_back(fragment->min_pos());
            seq2b[fragment->seq()].push_back(fragment->max_pos() + 1);
        }
    }
    BOOST_FOREACH (Seq2Boundaries::value_type& seq_and_boundaries, seq2b) {
        Boundaries& boundaries = seq_and_boundaries.second;
        select_boundaries(boundaries, min_distance());
    }
    bool result = false;
    BOOST_FOREACH (Block* block, *block_set()) {
        BOOST_FOREACH (Fragment* f, *block) {
            const Boundaries& boundaries = seq2b[f->seq()];
            size_t min_pos = nearest_element(boundaries, f->min_pos());
            BOOST_ASSERT(std::abs(int(min_pos - f->min_pos())) <
                         min_distance());
            size_t max_pos = nearest_element(boundaries, f->max_pos() + 1) - 1;
            BOOST_ASSERT(std::abs(int(max_pos - f->max_pos())) <
                         min_distance());
            if (min_pos != f->min_pos() || max_pos != f->max_pos()) {
                f->set_min_pos(min_pos);
                f->set_max_pos(max_pos);
                result = true;
            }
        }
    }
    return result;
}

const char* StickBoundaries::name_impl() const {
    return "Turn nearby fragment boundaries into one";
}

}

