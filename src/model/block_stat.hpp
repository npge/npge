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
class AlignmentStat {
public:
    /** Default constructor */
    AlignmentStat();

    /** Destructor */
    virtual ~AlignmentStat();

    /** Non-empty ident columns without gaps */
    int ident_nogap() const;

    /** Non-empty ident columns with gaps */
    int ident_gap() const;

    /** Non-empty non-ident columns without gaps */
    int noident_nogap() const;

    /** Non-empty non-ident columns with gaps */
    int noident_gap() const;

    /** Empty columns (consist only of gaps) */
    int pure_gap() const;

    /** All columns.
    Must be equal to sum of above variables.
    */
    int total() const;

    /** (max - min) / avg fragment length */
    float spreading() const;

    /** Number of rows with alignment */
    int alignment_rows() const;

    /** Minimal length of fragment */
    int min_fragment_length() const;

    /** Number of fragments overlapping their neighbours */
    int overlapping_fragments() const;

private:
    class Impl;

    Impl* impl_;

    friend void make_stat(AlignmentStat&, const Block*);
};

/** Make alignment stat of alignment */
void make_stat(AlignmentStat& stat, const Block* block);

/** Return proportion of columns, composed of size() equal letters.
If a fragment doesn't have alignment row attached,
then it is taken as is.

Gap columns are counted as half of non-gap columns:
(i_ng + 0.5 * i_g) / (i_ng + 0.5 * i_g + ni_ng + 0.5 * ni_g)

Column notation:
 - i = ident
 - g = gap
 - n = not
*/
float block_identity(const AlignmentStat& stat);

}

#endif

