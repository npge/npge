/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "StickBoundaries.hpp"
#include "Exception.hpp"
#include "stick_impl.hpp"
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
    if (vm.count("min-distance")) {
        int min_distance = vm["min-distance"].as<int>();
        if (min_distance < 0) {
            throw Exception("'min-distance' must be >= 0");
        }
        set_min_distance(min_distance);
    }
}

bool StickBoundaries::run_impl() const {
    Seq2Boundaries sb;
    bs_to_sb(sb, *block_set());
    stick_boundaries(sb, min_distance());
    return stick_fragments(*block_set(), sb, min_distance());
}

const char* StickBoundaries::name_impl() const {
    return "Turn nearby fragment boundaries into one";
}

}

