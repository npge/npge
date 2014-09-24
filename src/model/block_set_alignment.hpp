/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_BLOCK_SET_ALIGNMENT_HPP_
#define NPGE_BLOCK_SET_ALIGNMENT_HPP_

#include <iosfwd>
#include <map>
#include <vector>

#include "global.hpp"

namespace npge {

/** One row of blockset alignment */
struct BSRow {
    /** Constructor */
    BSRow();

    /** Ori of row in relation to sequence */
    int ori;

    /** Array of fragments, 0 for gap */
    Fragments fragments;
};

/** Blockset alignment */
typedef std::map<Sequence*, BSRow> BSA;

/** Vector of blockset alignment */
typedef std::vector<BSA> BSAs;

/** Return length of blockset alignment */
int bsa_length(const BSA& bsa);

/** Return is all sequences from bsa are circular */
bool bsa_is_circular(const BSA& bsa);

}

#endif

