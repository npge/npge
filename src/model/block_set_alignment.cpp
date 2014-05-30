/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "block_set_alignment.hpp"
#include "Sequence.hpp"
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

bool bsa_is_circular(const BSA& bsa) {
    BOOST_FOREACH (const BSA::value_type& seq_and_row, bsa) {
        Sequence* seq = seq_and_row.first;
        if (!seq->circular()) {
            return false;
        }
    }
    return true;
}

}

