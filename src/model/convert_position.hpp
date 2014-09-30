/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_CONVERT_POSITION_HPP_
#define NPGE_CONVERT_POSITION_HPP_

#include "global.hpp"

namespace npge {

/** Return block pos, corresponding to given fragment pos and block length */
pos_t block_pos(const Fragment* f, pos_t f_pos,
                pos_t block_length);

/** Return fragment pos, corresponding to given block pos and block length.
\see AlignmentRow::nearest_in_fragment
*/
pos_t fragment_pos(const Fragment* f, pos_t block_pos,
                   pos_t block_length);

/** Return sequence pos, corresponding to given fragment pos */
pos_t frag_to_seq(const Fragment* f, pos_t fragment_pos);

/** Return fragment pos, corresponding to given sequence pos */
int seq_to_frag(const Fragment* f, pos_t seq_pos);

/** Find columns at which slice is located */
void find_slice(pos_t& min_col, pos_t& max_col,
                const Block* host, const Block* slice);

}

#endif

