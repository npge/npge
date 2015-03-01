/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "convert_position.hpp"
#include "proportion.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "throw_assert.hpp"

namespace npge {

pos_t block_pos(const Fragment* f, pos_t f_pos,
                pos_t block_length) {
    ASSERT_GT(f->length(), 0);
    ASSERT_GT(block_length, 0);
    ASSERT_LTE(0, f_pos);
    ASSERT_LTE(f_pos, block_length);
    pos_t result;
    if (f->row()) {
        result = f->row()->map_to_alignment(f_pos);
        if (result == -1) {
            // FIXME ??
            if (f_pos < block_length / 2) {
                result = 0;
            } else {
                result = block_length;
            }
        }
    } else {
        result = proportion(f_pos, f->length(), block_length);
    }
    ASSERT_LTE(0, result);
    ASSERT_LTE(result, block_length);
    return result;
}

pos_t fragment_pos(const Fragment* f, pos_t block_pos,
                   pos_t block_length) {
    ASSERT_GT(f->length(), 0);
    ASSERT_GT(block_length, 0);
    ASSERT_LTE(0, block_pos);
    ASSERT_LTE(block_pos, block_length);
    pos_t result;
    if (f->row()) {
        result = f->row()->nearest_in_fragment(block_pos);
        if (result == -1) {
            // FIXME ??
            if (block_pos < block_length / 2) {
                result = 0;
            } else {
                result = f->length();
            }
        }
    } else {
        result = proportion(block_pos, block_length, f->length());
    }
    ASSERT_LTE(0, result);
    ASSERT_LTE(result, f->length());
    return result;
}

pos_t frag_to_seq(const Fragment* f, pos_t fragment_pos) {
    return f->begin_pos() + f->ori() * fragment_pos;
}

pos_t seq_to_frag(const Fragment* f, pos_t seq_pos) {
    return int(seq_pos - f->begin_pos()) * f->ori();
}

void find_slice(pos_t& min_col, pos_t& max_col,
                const Block* host, const Block* slice) {
    min_col = -1;
    max_col = -1;
    pos_t host_length = host->alignment_length();
    BOOST_FOREACH (const Fragment* s_f, *slice) {
        const Fragment* h_f = 0;
        BOOST_FOREACH (const Fragment* f, *host) {
            if (s_f->is_subfragment_of(*f)) {
                h_f = f;
                break;
            }
        }
        ASSERT_MSG(h_f, (host->name() + " " + slice->name() +
                         " " + s_f->id()).c_str());
        pos_t f1 = seq_to_frag(h_f, s_f->min_pos());
        pos_t f2 = seq_to_frag(h_f, s_f->max_pos());
        pos_t p1 = block_pos(h_f, f1, host_length);
        pos_t p2 = block_pos(h_f, f2, host_length);
        if (min_col == -1 || p1 < min_col) {
            min_col = p1;
        }
        if (min_col == -1 || p2 < min_col) {
            min_col = p2;
        }
        if (max_col == -1 || p1 > max_col) {
            max_col = p1;
        }
        if (max_col == -1 || p2 > max_col) {
            max_col = p2;
        }
    }
}

}

