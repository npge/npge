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

    /** Set all parameters to the most permissive */
    void allow_everything();

    /** Get min length of fragment */
    int min_fragment_length() const {
        return min_fragment_length_;
    }

    /** Set min length of fragment */
    void set_min_fragment_length(int min_fragment_length) {
        min_fragment_length_ = min_fragment_length;
    }

    /** Get max length of fragment (-1 = everything, default) */
    int max_fragment_length() const {
        return max_fragment_length_;
    }

    /** Set max length of fragment */
    void set_max_fragment_length(int max_fragment_length) {
        max_fragment_length_ = max_fragment_length;
    }

    /** Get min size of block */
    int min_block_size() const {
        return min_block_size_;
    }

    /** Set min size of block */
    void set_min_block_size(int min_block_size) {
        min_block_size_ = min_block_size;
    }

    /** Get max size of block (-1 = everything, default) */
    int max_block_size() const {
        return max_block_size_;
    }

    /** Set max size of block */
    void set_max_block_size(int max_block_size) {
        max_block_size_ = max_block_size;
    }

    /** Min fragment length spreading ((max - min) / avg).
    Defaults to 0 (everything).
    */
    float min_spreading() const {
        return min_spreading_;
    }

    /** Set min fragment length spreading.
    Defaults to 0.2.
    */
    void set_min_spreading(float min_spreading) {
        min_spreading_ = min_spreading;
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

    /** Return min block_identity().
    Should be applied only if alignment rows are known.
    */
    float min_identity() const {
        return min_identity_;
    }

    /** Set min block_identity().
    Defaults to 0.9.
    */
    void set_min_identity(float min_identity) {
        min_identity_ = min_identity;
    }

    /** Return max block_identity().
    Should be applied only if alignment rows are known.
    */
    float max_identity() const {
        return max_identity_;
    }

    /** Set max block_identity().
    Defaults to 1.0
    */
    void set_max_identity(float max_identity) {
        max_identity_ = max_identity;
    }

    /** Return min percentage of columns with gaps.
    Should be applied only if alignment rows are known.
    */
    float min_gaps() const {
        return min_gaps_;
    }

    /** Set min percentage of columns with gaps.
    Defaults to 0 (everything).
    */
    void set_min_gaps(float min_gaps) {
        min_gaps_ = min_gaps;
    }

    /** Return max percentage of columns with gaps.
    Should be applied only if alignment rows are known.
    */
    float max_gaps() const {
        return max_gaps_;
    }

    /** Set max percentage of columns with gaps.
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
    int max_fragment_length_;
    int min_block_size_;
    int max_block_size_;
    float min_spreading_; // fragment length spreading ((max - min) / avg)
    float max_spreading_; // fragment length spreading ((max - min) / avg)
    float min_identity_; // only if alignment rows are known
    float max_identity_; // only if alignment rows are known
    float min_gaps_; // only if alignment rows are known
    float max_gaps_; // only if alignment rows are known
};

}

#endif

