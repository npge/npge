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
    double spreading() const;

    /** Number of rows with alignment */
    int alignment_rows() const;

    /** Minimal length of fragment */
    int min_fragment_length() const;

    /** Number of fragments overlapping their neighbours */
    int overlapping_fragments() const;

    /** Number of corresponding nucleotides (a,t,g,c) */
    int letter_count(char letter) const;

    /** Part of G and C */
    double gc() const;

private:
    class Impl;

    Impl* impl_;

    friend void make_stat(AlignmentStat&, const Block*, int, int);
};

/** Make alignment stat of alignment.
\param stat Alignment stats
\param block Block
\param start first column to consider
\param stop last column to consider (-1 means last column of alignment)

start and stop affect gaps, ident and ATGC counters only.
*/
void make_stat(AlignmentStat& stat, const Block* block, int start = 0,
               int stop = -1);

/** Test one column of block.
\param block Block
\param column Index of column
\param ident Return if the column is identical
\param gap Return if the column contains gaps
\param pure_gap Return if the column contains only gaps
\param atgc Adds number of letters to the array.
    atgc is int[LETTERS_NUMBER]. Letters are converted into
    integers using char_to_size.

pure_gap columns are considered identical.
*/
void test_column(const Block* block, int column,
                 bool& ident, bool& gap, bool& pure_gap, int* atgc);

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
double block_identity(const AlignmentStat& stat);

/** Return proportion ident columns.
Gap columns are counted as half-of-column.
*/
double block_identity(int ident_nogap, int ident_gap,
                      int noident_nogap, int noident_gap);

}

#endif

