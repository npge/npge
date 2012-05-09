/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PAIR_ALIGNER_HPP
#define BR_PAIR_ALIGNER_HPP

#include <vector>

#include "global.hpp"

namespace bloomrepeats {

class PairAligner {
public:
    /** Constructor.
    \param max_errors Max number of errors in pair alignment,
        that may happen on a batch.
    \param gap_range Max distance from main diagonal of considered
        states of pair alignment. The more gap_range, the more time.
    */
    PairAligner(int max_errors = 10, int gap_range = 5);

    void set_first(const char* start, int size);

    void set_second(const char* start, int size);

    int gap_range() const {
        return gap_range_;
    }

    void set_gap_range(int gap_range);

    int max_errors() const {
        return max_errors_;
    }

    void set_max_errors(int max_errors);

    bool no_tail() const {
        return no_tail_;
    }

    void set_no_tail(bool no_tail) {
        no_tail_ = no_tail;
    }

    void align(int& first_last, int& second_last) const;

private:
    mutable std::vector<int> matrix_;
    int gap_range_, max_errors_;
    const char* first_start_;
    const char* second_start_;
    int first_size_, second_size_;
    bool no_tail_;

    void cut_tail(int& first_last, int& second_last) const;

    int rows() const {
        return first_size_;
    }

    int cols() const {
        return second_size_;
    }

    int side() const;

    int row_size() const;

    int min_col(int row) const;

    int max_col(int row) const;

    int& at(int row, int col) const;

    bool in(int row, int col) const;

    int substitution(int row, int col) const;

    void adjust_matrix_size();
};

}

#endif

