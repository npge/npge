/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_HIT_HPP_
#define BR_HIT_HPP_

#include "FragmentCollection.hpp"

namespace bloomrepeats {

/** Set of fragments */
typedef std::set<Fragment*, FragmentCompare> FragmentsSet;

/** Collection of fragments */
typedef FragmentCollection<Fragment*, FragmentsSet> S2F;

/** Return if hit is internal for block set (s2f).
This means, each fragment of hit is a subfragment of one fragment
from s2f and size of corresponding block from blockset is
less than size of hit.
*/
bool is_internal_hit(const S2F& s2f, const Block* hit);

/** Return if block contain self-overlapping fragments */
bool has_self_overlaps(Block* block);

/** Shorten or remove self-overlapping fragments */
void fix_self_overlaps(Block* block);

}

#endif
