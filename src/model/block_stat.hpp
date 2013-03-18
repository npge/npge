/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_STAT_HPP_
#define BR_BLOCK_STAT_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Column stat of alignment */
struct AlignmentStat {
    /** Default constructor */
    AlignmentStat();

    /** Non-empty ident columns without gaps */
    int ident_nogap;

    /** Non-empty ident columns with gaps */
    int ident_gap;

    /** Non-empty non-ident columns without gaps */
    int noident_nogap;

    /** Non-empty non-ident columns with gaps */
    int noident_gap;

    /** Empty columns (consist only of gaps) */
    int pure_gap;

    /** All columns.
    Must be equal to sum of above variables.
    */
    int total;

    /** (max - min) / avg fragment length */
    float spreading;
};

/** Make alignment stat of alignment */
void make_stat(AlignmentStat& stat, const Block* block);

/** Return proportion of columns, composed of size() equal letters.
If a fragment doesn't have alignment row attached,
then it is taken as is.
*/
float block_identity(const AlignmentStat& stat, bool allow_gaps = false);

}

#endif

