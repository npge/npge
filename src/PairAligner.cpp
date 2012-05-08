/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/assert.hpp>

#include "PairAligner.hpp"

namespace bloomrepeats {

PairAligner::PairAligner(int max_errors, int gap_range):
    gap_range_(gap_range), max_errors_(max_errors),
    first_start_(0), second_start_(0),
    first_size_(0), second_size_(0),
    no_tail_(true) {
    gap_range_ = std::min(gap_range_, max_errors_);
}

void PairAligner::set_first(const char* start, int size) {
    first_start_ = start;
    first_size_ = size;
    adjust_matrix_size();
}

void PairAligner::set_second(const char* start, int size) {
    second_start_ = start;
    second_size_ = size;
    adjust_matrix_size();
}

void PairAligner::set_gap_range(int gap_range) {
    gap_range_ = gap_range;
    max_errors_ = std::max(gap_range_, max_errors_);
    adjust_matrix_size();
}

void PairAligner::set_max_errors(int max_errors) {
    max_errors_ = max_errors;
    gap_range_ = std::min(gap_range_, max_errors_);
    adjust_matrix_size();
}

void PairAligner::align(int& r_row, int& r_col) const {
    at(0, 0) = 0;
    for (int row = 1; row <= gap_range_; row++) {
        at(row, 0) = row;
    }
    for (int col = 1; col <= gap_range_; col++) {
        at(0, col) = col;
    }
    for (int row = 1; row < side(); row++) {
        int start_col = std::max(1, min_col(row));
        int min_score_col = start_col;
        for (int col = start_col; col <= max_col(row); col++) {
            BOOST_ASSERT(col >= 1 && col < side());
            BOOST_ASSERT(in(row, col));
            int score = at(row - 1, col - 1) + substitution(row, col);
            if (in(row - 1, col)) {
                score = std::min(score, at(row - 1, col) + 1);
            }
            if (in(row, col - 1)) {
                score = std::min(score, at(row, col - 1) + 1);
            }
            at(row, col) = score;
            if (score < at(row, min_score_col)) {
                min_score_col = col;
            }
        }
        if (at(row, min_score_col) > max_errors_) {
            r_row = row;
            r_col = min_score_col;
            break;
        }
    }
    r_col = std::min(cols() - 1, r_col);
    r_row = std::min(rows() - 1, r_row);
    if (no_tail_) {
        cut_tail(r_row, r_col);
    }
}

void PairAligner::cut_tail(int& r_row, int& r_col) const {
    while (true) {
        if (in(r_row - 1, r_col) && at(r_row - 1, r_col) < at(r_row, r_col)) {
            r_row -= 1;
        } else if (in(r_row, r_col - 1) &&
                   at(r_row, r_col - 1) < at(r_row, r_col)) {
            r_col -= 1;
        } else if (in(r_row - 1, r_col - 1) &&
                   at(r_row - 1, r_col - 1) < at(r_row, r_col)) {
            r_row -= 1;
            r_col -= 1;
        } else {
            break;
        }
    }
}

int PairAligner::side() const {
    return std::min(rows(), cols()) + gap_range_;
}

int PairAligner::row_size() const {
    return 1 + 2 * gap_range_;
}

int PairAligner::min_col(int row) const {
    return std::max(0, row - gap_range_);
}

int PairAligner::max_col(int row) const {
    return std::min(side() - 1, row + gap_range_);
}

int& PairAligner::at(int row, int col) const {
    return matrix_[row * row_size() + gap_range_ + col - row];
}

bool PairAligner::in(int row, int col) const {
    return row >= 0 && row < side() &&
           col >= min_col(row) && col <= max_col(row);
}

int PairAligner::substitution(int row, int col) const {
    return first_start_[row] == second_start_[col] ? 0 : 1;
}

void PairAligner::adjust_matrix_size() {
    matrix_.resize(side() * row_size());
}

}

