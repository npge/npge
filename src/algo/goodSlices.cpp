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
    std::vector<int> score_sum_; // prefix sum
    int min_length_;
    int min_end_;
    int min_ident_;
    int block_length_;

public:
    GoodSlicer(const Scores& score,
               int min_length, int min_end, int min_ident):
        min_length_(min_length),
        min_end_(min_end),
        min_ident_(min_ident) {
        block_length_ = score.size();
        score_sum_.resize(block_length_ + 1);
        score_sum_[0] = 0;
        for (int i = 0; i < block_length_; i++) {
            score_sum_[i + 1] = score_sum_[i] + score[i];
        }
    }

    int countScore(int start, int stop) const {
        return score_sum_[stop + 1] - score_sum_[start - 1 + 1];
    }

    bool goodSlice(int start) const {
        int stop = start + min_length_ - 1;
        return countScore(start, stop) >= min_ident_;
    }

    bool goodLeftEnd(int start) const {
        int stop = start + min_end_ - 1;
        return countScore(start, stop) >= min_end_;
    }

    bool goodRightEnd(int stop) const {
        int start = stop - min_end_ + 1;
        return countScore(start, stop) >= min_end_;
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
                start1 + min_end_ - 1 < stop1) {
            start1 = start1 + 1;
        }
        while (!goodRightEnd(stop1) &&
                start1 + min_end_ - 1 < stop1) {
            stop1 = stop1 - 1;
        }
        return StartStop(start1, stop1);
    }

    bool valid(const StartStop& self) const {
        return ssLength(self) >= min_length_ &&
            self.first >= 0 && self.second < block_length_;
    }

    // Return list of joined slices
    Coordinates joinedSlices() const {
        Coordinates slices0;
        bool prev_good = false;
        for (int i = 0; i <= block_length_ - min_length_; i++) {
            bool curr_good = goodSlice(i);
            if (curr_good) {
                if (prev_good) {
                    // increase previous slice
                    ASSERT_GT(slices0.size(), 0);
                    slices0.back().second += 1;
                } else {
                    // add new slice
                    int stop = i + min_length_ - 1;
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
        if (min_length_ > block_length_ || min_length_ <= 0) {
            return false;
        }
        if (min_end_ > min_length_ || min_end_ < 0) {
            return false;
        }
        if (min_ident_ > min_length_ || min_ident_ < 0) {
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
            if (ssLength(selected) >= min_length_) {
                result.push_back(selected);
                slices = excludeSlice(slices, selected);
            }
        }
        return result;
    }
};

Coordinates goodSlices(const Scores& score, int min_length,
                       int min_end, int min_ident) {
    GoodSlicer slicer(score, min_length,
                      min_end, min_ident);
    return slicer.calculate();
}

}
