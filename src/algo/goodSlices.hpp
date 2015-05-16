/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GOOD_SLICES_HPP_
#define NPGE_GOOD_SLICES_HPP_

#include <vector>

namespace npge {

typedef std::pair<int, int> StartStop; // start, stop
typedef std::vector<StartStop> Coordinates;
typedef std::vector<bool> Columns;

Coordinates goodSlices(const Columns& columns, int min_length,
                       int min_end, int min_ident);

}

#endif
