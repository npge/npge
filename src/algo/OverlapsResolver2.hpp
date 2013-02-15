/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OVERLAPS_RESOLVER2_HPP_
#define BR_OVERLAPS_RESOLVER2_HPP_

#include "OverlapsResolver.hpp"

namespace bloomrepeats {

/** Resolve overlaping fragments (version 2).
\warning Source blocks are taken from "other" block set.
*/
class OverlapsResolver2 : public OverlapsResolver {
public:
    /** Constructor */
    OverlapsResolver2(int min_distance = 30);

    /** Get min distance between fragment boundaries */
    int min_distance() const {
        return min_distance_;
    }

    /** Set min distance between fragment boundaries */
    void set_min_distance(int min_distance) {
        min_distance_ = min_distance;
    }

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;

private:
    int min_distance_;
};

}

#endif

