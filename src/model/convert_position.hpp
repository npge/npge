/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONVERT_POSITION_HPP_
#define BR_CONVERT_POSITION_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Return block pos, corresponding to given fragment pos and block length */
int block_pos(const Fragment* f, int f_pos, int block_length);

/** Return fragment pos, corresponding to given block pos and block length.
\see AlignmentRow::nearest_in_fragment
*/
int fragment_pos(const Fragment* f, int block_pos, int block_length);

/** Return sequence pos, corresponding to given fragment pos */
size_t frag_to_seq(const Fragment* f, int fragment_pos);

/** Return fragment pos, corresponding to given sequence pos */
int seq_to_frag(const Fragment* f, size_t seq_pos);

}

#endif

