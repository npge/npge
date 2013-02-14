/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "StickBoundaries.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Exception.hpp"
#include "boundaries.hpp"
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

