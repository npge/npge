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

StickBoundaries::StickBoundaries(int min_distance) {
    add_opt("min-distance",
            "Min distance between fragment boundaries",
            min_distance);
    add_opt_rule("min-distance >= 0");
    declare_bs("target", "Target blockset");
}

void StickBoundaries::run_impl() const {
    int min_distance = opt_value("min-distance").as<int>();
    Seq2Boundaries sb;
    bs_to_sb(sb, *block_set());
    stick_boundaries(sb, min_distance);
    stick_fragments(*block_set(), sb, min_distance);
}

const char* StickBoundaries::name_impl() const {
    return "Turn nearby fragment boundaries into one";
}

}

