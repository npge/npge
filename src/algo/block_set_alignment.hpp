/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_SET_ALIGNMENT_HPP_
#define BR_BLOCK_SET_ALIGNMENT_HPP_

#include <map>
#include <vector>

#include "global.hpp"

namespace bloomrepeats {

/** Array of fragments, 0 for gap */
typedef std::vector<Fragment*> Fragments;

/** One row of block set alignment */
struct BSRow {
    /** Constructor */
    BSRow();

    /** Ori of row in relation to sequence */
    int ori;

    /** Array of fragments, 0 for gap */
    Fragments fragments;
    // TODO circular genomes
    // TODO start of row in arbitrary part of sequence
};

/** Block set alignment */
typedef std::map<Sequence*, BSRow> BSA;

/** Return length of block set alignment */
int bsa_length(const BSA& bsa);

/** Create blocks set alignment row of the sequence.
Output BSA is not an alignment.
If rows is empty, all sequences of the chromosome will be added.
If rows is not empty, only sequences existing in rows will be added.
*/
void bsa_make_rows(BSA& rows, const BlockSet& bs,
                   const std::string& chromosome);

/** Inverse alignment */
void bsa_inverse(BSA& aln);

/** Create blocks set alignment row of the sequence */
void bsa_align(BSA& both, int& score,
               const BSA& first, const BSA& second);

/** Produce alignment from map of trivial rows */
void bsa_make_aln(BSA& aln, const BSA& rows);
// TODO use tree

}

#endif

