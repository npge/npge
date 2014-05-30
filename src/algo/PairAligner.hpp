/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PAIR_ALIGNER_HPP
#define BR_PAIR_ALIGNER_HPP

#include <vector>
#include <string>

#include "global.hpp"

namespace npge {

/** Find the end of good alignment using Needleman-Wunsch with gap frame */
class PairAligner {
public:
    /** Constructor.
    \param max_errors Max number of errors in pair alignment.
        Alignment stops, when max errors accumulated.
    \param gap_range Max distance from main diagonal of considered
        states of pair alignment. The more gap_range, the more time.
    \param gap_penalty gap open or extension penalty.
    */
    PairAligner();

    /** Destructor */
    virtual ~PairAligner();

    /** Return a pointer to global thread-local default PairAligner */
    static PairAligner* default_aligner();

    /** Set first sequence */
    void set_first(const char* start, int size);

    /** Set second sequence */
    void set_second(const char* start, int size);

    /** Get gap range */
    int gap_range() const;

    /** Set gap range.
    \code gap frame = 2 * gap range + 1 \endcode
    Gap range is max distance from main diagonal of considered
    states of pair alignment. The more gap_range, the more time.
    */
    void set_gap_range(int gap_range);

    /** Get max errors */
    int max_errors() const;

    /** Set max errors.
    Mismatch or gap are considered as 1 error.
    Alignment stops, when max errors accumulated.
    */
    void set_max_errors(int max_errors);

    /** Get gap penalty */
    int gap_penalty() const;

    /** Set gap penalty */
    void set_gap_penalty(int gap_penalty);

    /** Get mismatch penalty */
    int mismatch_penalty() const;

    /** Set mismatch penalty */
    void set_mismatch_penalty(int mismatch_penalty);

    /** Get whether bad alignment tail would be stripped out */
    bool no_tail() const;

    /** Set whether bad alignment tail would be stripped out.
    Consider the following alignment:
    \verbatim
    TTCCGGTGCTGCGaggga
    TTCCGGTGCTGCGcctct
    \endverbatim
    After tail stripping, it would be changed to
    \verbatim
    TTCCGGTGCTGCG
    TTCCGGTGCTGCG
    \endverbatim

    Defaults to \c true.
    */
    void set_no_tail(bool no_tail);

    /** Return if the alignment is local */
    bool local() const;

    /** Set if the alignment is local.
    If local alignment is used, then max_errors() must be -1.
    Default: false.
    */
    void set_local(bool local);

    /** Run alignment algorithm.
    \param first_last Last aligned position in first sequence (output)
    \param second_last Last aligned position in second sequence (output)
    \param first_str Alignment row (like "a-tgc") (optional output)
    \param second_str Alignment row (like "a-tgc") (optional output)
    \param alignment Array of pairs of positions, -1 as gaps (optional output)
    \param gap Char used as gap for first_str and second_str
    */
    void align(int& first_last, int& second_last,
               std::string* first_str = 0, std::string* second_str = 0,
               PairAlignment* alignment = 0,
               char gap = '-') const;

    /** Return if two sequences can be globally aligned.
    If both \p first_last and \p second_last are provided,
    the result of align() is written to them.

    \note This method calls set_first() and set_second().
    */
    bool aligned(const std::string& first, const std::string& second,
                 int* first_last = 0, int* second_last = 0);

private:
    class Impl;
    Impl* impl_;

    void export_alignment(int row, int col,
                          std::string* first_str,
                          std::string* second_str,
                          PairAlignment* alignment,
                          char gap) const;

    int substitution(int row, int col) const;
};

}

#endif

