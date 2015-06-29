/* Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "goodSlices.hpp"
#include "throw_assert.hpp"

namespace npge {

static int ssLength(const StartStop& ss) {
    return ss.second - ss.first + 1;
}

class GoodSlicer {
private:
    std::vector<int> score_;
    std::vector<int> score_sum_; // prefix sum
    int frame_length_;
    int end_length_;
    int frame_score_;
    int end_score_;
    int block_length_;

public:
    GoodSlicer(const Scores& score,
               int frame_length, int end_length,
               int min_identity):
        score_(score),
        frame_length_(frame_length),
        end_length_(end_length),
        frame_score_(frame_length * min_identity),
        end_score_(end_length * min_identity) {
        block_length_ = score.size();
        score_sum_.resize(block_length_ + 1);
        score_sum_[0] = 0;
        for (int i = 0; i < block_length_; i++) {
            int value = score[i];
            // for gap columns, multiply score by min_ident
            value = (value == MAX_COLUMN_SCORE) ? value :
                (value * min_identity / MAX_COLUMN_SCORE);
            score_sum_[i + 1] = score_sum_[i] + value;
        }
    }

    int countScore(int start, int stop) const {
        return score_sum_[stop + 1] - score_sum_[start - 1 + 1];
    }

    bool goodSlice(int start) const {
        int stop = start + frame_length_ - 1;
        return countScore(start, stop) >= frame_score_;
    }

    bool goodLeftEnd(int start) const {
        int stop = start + end_length_ - 1;
        return score_[start] == MAX_COLUMN_SCORE &&
               countScore(start, stop) >= end_score_;
    }

    bool goodRightEnd(int stop) const {
        int start = stop - end_length_ + 1;
        return score_[stop] == MAX_COLUMN_SCORE &&
               countScore(start, stop) >= end_score_;
    }

    bool overlaps(const StartStop& self,
                  const StartStop& other) const {
        if (other.first <= self.first &&
                self.first <= other.second) {
            return true;
        }
        if (self.first <= other.first &&
                other.first <= self.second) {
            return true;
        }
        return false;
    }

    StartStop exclude(const StartStop& self,
                      const StartStop& other) const {
        int start1 = self.first;
        int stop1 = self.second;
        if (other.first <= self.first &&
                self.first <= other.second) {
            start1 = other.second + 1;
        }
        if (other.first <= self.second &&
                self.second <= other.second) {
            stop1 = other.first - 1;
        }
        return StartStop(start1, stop1);
    }

    StartStop strip(const StartStop& self) const {
        if (!valid(self)) {
            return self;
        }
        int start1 = self.first;
        int stop1 = self.second;
        while (!goodLeftEnd(start1) &&
                start1 + end_length_ - 1 < stop1) {
            start1 = start1 + 1;
        }
        while (!goodRightEnd(stop1) &&
                start1 + end_length_ - 1 < stop1) {
            stop1 = stop1 - 1;
        }
        return StartStop(start1, stop1);
    }

    bool valid(const StartStop& self) const {
        return ssLength(self) >= frame_length_ &&
            self.first >= 0 && self.second < block_length_;
    }

    // Return list of joined slices
    Coordinates joinedSlices() const {
        Coordinates slices0;
        bool prev_good = false;
        int max_i = block_length_ - frame_length_;
        for (int i = 0; i <= max_i; i++) {
            bool curr_good = goodSlice(i);
            if (curr_good) {
                if (prev_good) {
                    // increase previous slice
                    ASSERT_GT(slices0.size(), 0);
                    slices0.back().second += 1;
                } else {
                    // add new slice
                    int stop = i + frame_length_ - 1;
                    slices0.push_back(StartStop(i, stop));
                }
            }
            prev_good = curr_good;
        }
        Coordinates slices;
        BOOST_FOREACH (const StartStop& slice, slices0) {
            StartStop slice1 = strip(slice);
            if (valid(slice1)) {
                slices.push_back(slice1);
            }
        }
        return slices;
    }

    StartStop maxSlice(const Coordinates& slices) const {
        StartStop result = slices.front();
        BOOST_FOREACH (const StartStop& slice, slices) {
            if (ssLength(slice) > ssLength(result)) {
                result = slice;
            }
        }
        return result;
    }

    // Exclude selected slice from slices, return new slices
    Coordinates excludeSlice(const Coordinates& slices,
                             const StartStop& selected) const {
        Coordinates slices1;
        BOOST_FOREACH (const StartStop& slice, slices) {
            if (!overlaps(slice, selected)) {
                slices1.push_back(slice);
            } else {
                StartStop slice1 = exclude(slice, selected);
                slice1 = strip(slice1);
                if (valid(slice1)) {
                    slices1.push_back(slice1);
                }
            }
        }
        return slices1;
    }

    bool parametersAreCorrect() const {
        if (frame_length_ > block_length_ ||
                frame_length_ <= 0) {
            return false;
        }
        if (end_length_ > frame_length_ || end_length_ < 0) {
            return false;
        }
        return true;
    }

    Coordinates calculate() const {
        if (!parametersAreCorrect()) {
            return Coordinates();
        }
        Coordinates slices = joinedSlices();
        Coordinates result;
        while (!slices.empty()) {
            StartStop selected = maxSlice(slices);
            if (ssLength(selected) >= frame_length_) {
                result.push_back(selected);
                slices = excludeSlice(slices, selected);
            }
        }
        return result;
    }
};

Coordinates goodSlices(const Scores& score,
                       int frame_length, int end_length,
                       int min_identity) {
    GoodSlicer slicer(score,
                      frame_length, end_length,
                      min_identity);
    return slicer.calculate();
}

}
