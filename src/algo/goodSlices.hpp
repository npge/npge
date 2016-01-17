/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GOOD_SLICES_HPP_
#define NPGE_GOOD_SLICES_HPP_

#include <vector>

#include "goodColumns.hpp"

namespace npge {

typedef std::pair<int, int> StartStop; // start, stop
typedef std::vector<StartStop> Coordinates;

Coordinates goodSlices(const Scores& score,
                       int frame_length, int end_length,
                       int min_identity, int min_length);

}

#endif
