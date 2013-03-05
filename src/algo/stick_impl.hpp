/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_STICK_IMPL_HPP_
#define BR_STICK_IMPL_HPP_

#include <map>

#include "global.hpp"
#include "SortedVector.hpp"
#include "Graph.hpp"
#include "boundaries.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

typedef std::pair<Sequence*, size_t> Point;
typedef SortedVector<Point> Points;
typedef std::pair<Point, Point> PointsPair;
typedef Graph<Point> PointsGraph;
typedef std::map<Sequence*, Boundaries> Seq2Boundaries;

/** Add all fragment boundaries to sb.
\note Boundaries are not sorted.
*/
void bs_to_sb(Seq2Boundaries& sb, const BlockSet& bs);

/** Return if two boundaries collections are equal */
bool sb_equal(const Seq2Boundaries& x, const Seq2Boundaries& y);

/** Return boundaries collection matches block set (same boundary sets) */
bool sb_match_bs(const Seq2Boundaries& sb, const BlockSet& bs);

/** Sort the vector and merge too close boundaries together */
void stick_boundaries(Seq2Boundaries& sb, int min_distance);

/** Change fragments of block set to boundaries from sb.
min_distance is used only for assertations.
*/
bool stick_fragments(BlockSet& bs, const Seq2Boundaries& sb, int min_distance);

}

#endif

