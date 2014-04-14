/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
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

namespace bloomrepeats {

int block_pos(const Fragment* f, int f_pos, int block_length) {
    ASSERT_GT(f->length(), 0);
    ASSERT_GT(block_length, 0);
    ASSERT_LTE(0, f_pos);
    ASSERT_LTE(f_pos, block_length);
    int result;
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

int fragment_pos(const Fragment* f, int block_pos, int block_length) {
    ASSERT_GT(f->length(), 0);
    ASSERT_GT(block_length, 0);
    ASSERT_LTE(0, block_pos);
    ASSERT_LTE(block_pos, block_length);
    int result;
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

size_t frag_to_seq(const Fragment* f, int fragment_pos) {
    return f->begin_pos() + f->ori() * fragment_pos;
}

int seq_to_frag(const Fragment* f, size_t seq_pos) {
    return int(seq_pos - f->begin_pos()) * f->ori();
}

void find_slice(int& min_col, int& max_col,
                const Block* host, const Block* slice) {
    min_col = -1;
    max_col = -1;
    int host_length = host->alignment_length();
    BOOST_FOREACH (const Fragment* s_f, *slice) {
        const Fragment* h_f = 0;
        BOOST_FOREACH (const Fragment* f, *host) {
            if (s_f->is_subfragment_of(*f)) {
                h_f = f;
                break;
            }
        }
        BOOST_ASSERT_MSG(h_f, (host->name() + " " + slice->name() +
                               " " + s_f->id()).c_str());
        int f1 = seq_to_frag(h_f, s_f->min_pos());
        int f2 = seq_to_frag(h_f, s_f->max_pos());
        int p1 = block_pos(h_f, f1, host_length);
        int p2 = block_pos(h_f, f2, host_length);
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

