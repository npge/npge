/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BOUNDARIES_HPP_
#define BR_BOUNDARIES_HPP_

#include "SortedVector.hpp"

namespace bloomrepeats {

/** Vector of positions.
Elements of the vector may be Fragment::min_pos() or Fragment::max_pos() + 1.
*/
typedef SortedVector<size_t> Boundaries;

/** Vector of floats */
typedef std::vector<float> Floats;

/** Return average value of the vector */
size_t avg_element(const Boundaries& boundaries);

/** Return average value of the vector */
double avg_element_double(const Boundaries& boundaries);

/** Return average value of the vector */
double avg_element_double(const Floats& floats);

/** Return value of the nearest element to the position.
The vector must be sorted in ascending.
*/
size_t nearest_element(const Boundaries& boundaries, size_t pos);

/** Sort the vector and merge too close elements together.
\param boundaries List of boundaries.
\param min_distance min distance between input boundaries,
    which guarantees that they will not be merged into one boundary.
\param length Length of sequence (boundaries are coordinates in the sequence).

If distance between a boundary and first/last nucleotide
is less than min_distance, then the boundary is moved to first/last nucleotide.
*/
void select_boundaries(Boundaries& boundaries, int min_distance, size_t length);

}

#endif

