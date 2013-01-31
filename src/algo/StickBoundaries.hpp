/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_STICK_BOUNDARIES_HPP_
#define BR_STICK_BOUNDARIES_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Move boundaries of fragments to turn nearby boundaries into one.
It is recommended to \ref Filter "filter" fragments before
applying this processor.

\warning Connections between fragments may be broken.
    They must be \ref Connector "re-connected" after applying this processor.

\warning Alignment of blocks may be broken.
*/
class StickBoundaries : public Processor {
public:
    /** Constructor */
    StickBoundaries(int min_distance = 100);

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

