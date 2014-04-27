/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FIND_LOW_SIMILAR_HPP_
#define BR_FIND_LOW_SIMILAR_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Find regions of low similarity in blocks */
class FindLowSimilar : public BlocksJobs {
public:
    /** Constructor */
    FindLowSimilar();

    /** Region of same color */
    struct Region {
        int start_;
        int stop_;
        int good_;
        int weight_;

        int length() const;
        void set_weight(int weight_factor);
    };

    /** List of regions (ordered by start) */
    typedef std::vector<Region> Regions;

    /** Return weight factor corresponding to min_identity */
    static int get_weight_factor(double min_identity);

    /** Create list of regions from list of color.
    Weight factor is applied to bad regions.
    */
    static Regions make_regions(const std::vector<bool>& good_col,
                                int weight_factor);

    /** Find index of shortest region (by weight) */
    static int find_min_region(const Regions& regions);

    /** Replace a region with color of its neighbors */
    static Regions merge_region(Regions& regions, int index);

    /** Reduce shortest region while its length < min_length */
    static void reduce_regions(Regions& regions, int min_length);

protected:
    ThreadData* before_thread_impl() const;
    void process_block_impl(Block* block, ThreadData*) const;
    void after_thread_impl(ThreadData* data) const;
    const char* name_impl() const;
};

}

#endif

