/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "convert_position.hpp"
#include "proportion.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"

namespace bloomrepeats {

int block_pos(const Fragment* f, int f_pos, int block_length) {
    BOOST_ASSERT(f->length() > 0);
    BOOST_ASSERT(block_length > 0);
    BOOST_ASSERT(0 <= f_pos && f_pos <= block_length);
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
    BOOST_ASSERT(0 <= result && result <= block_length);
    return result;
}

int fragment_pos(const Fragment* f, int block_pos, int block_length) {
    BOOST_ASSERT(f->length() > 0);
    BOOST_ASSERT(block_length > 0);
    BOOST_ASSERT(0 <= block_pos && block_pos <= block_length);
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
    BOOST_ASSERT(0 <= result && result <= f->length());
    return result;
}

size_t frag_to_seq(const Fragment* f, int fragment_pos) {
    return f->begin_pos() + f->ori() * fragment_pos;
}

int seq_to_frag(const Fragment* f, size_t seq_pos) {
    return int(seq_pos - f->begin_pos()) * f->ori();
}

}

