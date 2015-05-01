/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "block_set_alignment.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "throw_assert.hpp"

namespace npge {

BSRow::BSRow():
    ori(1) {
}

int bsa_length(const BSA& bsa) {
    if (bsa.empty()) {
        return 0;
    } else {
        int length = bsa.begin()->second.fragments.size();
        BOOST_FOREACH (const BSA::value_type& seq_and_row, bsa) {
            const BSRow& row = seq_and_row.second;
            ASSERT_EQ(row.fragments.size(), length);
        }
        return length;
    }
}

static Fragment* firstNonGap(const Fragments& ff) {
    BOOST_FOREACH (Fragment* f, ff) {
        if (f) {
            return f;
        }
    }
    return 0;
}

static Fragment* lastNonGap(const Fragments& ff) {
    BOOST_REVERSE_FOREACH (Fragment* f, ff) {
        if (f) {
            return f;
        }
    }
    return 0;
}

bool bsa_is_circular(const BSA& bsa) {
    BOOST_FOREACH (const BSA::value_type& seq_and_row, bsa) {
        Sequence* seq = seq_and_row.first;
        if (!seq->circular()) {
            return false;
        }
        const Fragments& ff = seq_and_row.second.fragments;
        Fragment* first = firstNonGap(ff);
        Fragment* last = lastNonGap(ff);
        if (!first || first->min_pos() != 0) {
            return false;
        }
        if (!last || last->max_pos() != seq->size() - 1) {
            return false;
        }
    }
    return true;
}

}

