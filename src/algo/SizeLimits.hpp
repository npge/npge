/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SIZE_LIMITS_HPP_
#define BR_SIZE_LIMITS_HPP_

#include "po.hpp"

namespace bloomrepeats {

/** Size limits for block and fragment */
class SizeLimits {
public:
    /** Constructor */
    SizeLimits(int min_fragment_length = 100, int min_block_size = 2);

    /** Get min length of fragment */
    int min_fragment_length() const {
        return min_fragment_length_;
    }

    /** Set min length of fragment */
    void set_min_fragment_length(int min_fragment_length) {
        min_fragment_length_ = min_fragment_length;
    }

    /** Get min size of block */
    int min_block_size() const {
        return min_block_size_;
    }

    /** Set min size of block */
    void set_min_block_size(int min_block_size) {
        min_block_size_ = min_block_size;
    }

    /** Max fragment length spreading ((max - min) / avg) */
    float max_spreading() const {
        return max_spreading_;
    }

    /** Set max fragment length spreading.
    Defaults to 0.2.
    */
    void set_max_spreading(float max_spreading) {
        max_spreading_ = max_spreading;
    }

    /** Return min identity on non-gap columns.
    Should be applied only if alignment rows are known.
    */
    float min_identity() const {
        return min_identity_;
    }

    /** Set min identity on non-gap columns.
    Defaults to 0.9.
    */
    void set_min_identity(float min_identity) {
        min_identity_ = min_identity;
    }

    /** Return max percentage of columns with gaps.
    Should be applied only if alignment rows are known.
    */
    float max_gaps() const {
        return max_gaps_;
    }

    /** Set min identity.
    Defaults to 0.2.
    */
    void set_max_gaps(float max_gaps) {
        max_gaps_ = max_gaps;
    }

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

private:
    int min_fragment_length_;
    int min_block_size_;
    float max_spreading_; // fragment length spreading ((max - min) / avg)
    float min_identity_; // only if alignment rows are known
    float max_gaps_; // only if alignment rows are known
};

}

#endif

