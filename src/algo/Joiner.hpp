/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_JOIN_APPROVER_HPP_
#define BR_JOIN_APPROVER_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace bloomrepeats {

/** Utility object, making decision whether blocks/fragments can be merged.
Blocks/fragments must be joinable (Block::can_join and Fragment::can_join).

\ref Block::weak() "Weak" blocks can't be joined.
*/
class Joiner : public Processor {
public:
    /** Constructor.
    \param max_dist Max allowed \ref Fragment::dist_to "distance".
        Value -1 means that this limitation is not applied.
    \param ratio_to_fragment Max allowed gap length to fragment length ratio.
        A negative number means that this limitation is not applied.
    \param gap_ratio Max allowed ratio of gaps' lengths (inside a block).
        A negative number means that this limitation is not applied.
    */
    Joiner(int max_dist = -1, float ratio_to_fragment = -1,
           float gap_ratio = -1);

    /** Get max allowed distance */
    int max_dist() const {
        return max_dist_;
    }

    /** Set max allowed distance */
    void set_max_dist(int max_dist) {
        max_dist_ = max_dist;
    }

    /** Get max allowed gap length to fragment length ratio */
    float ratio_to_fragment() const {
        return ratio_to_fragment_;
    }

    /** Set max allowed gap length to fragment length ratio */
    void set_ratio_to_fragment(float ratio_to_fragment) {
        ratio_to_fragment_ = ratio_to_fragment;
    }

    /** Get max allowed ratio of gaps' lengths (inside a block) */
    float gap_ratio() const {
        return gap_ratio_;
    }

    /** Set max allowed ratio of gaps' lengths (inside a block) */
    void set_gap_ratio(float gap_ratio) {
        gap_ratio_ = gap_ratio;
    }

    /** Return if these fragments can be joined (simple check).
    Fragments can be joined if they share the same sequence and ori
    and \ref is_neighbor "are neighbors".
    */
    static bool can_join(Fragment* one, Fragment* another);

    /** Return if these blocks can be joined (simple check).
     - 1, if fragments of 'one' should preceed fragments from 'another';
     - -1, if fragments of 'another' should preceed fragments from 'one';
     - 0, if blocks can't be joined.

    Blocks can be joined, if they match and all the fragments from
    the first block has an unique neighbor with the same ori
    from the second block.
    */
    static int can_join(Block* one, Block* another);

    /** Merge fragments and return new larger fragment */
    static Fragment* join(Fragment* one, Fragment* another);

    /** Return joined blocks, if these two blocks can be joined */
    static Block* join(Block* one, Block* another, int logical_ori);

    /** Return if two fragments can be joined.
    The fragments must be joinable (Fragment::can_join) and non empty.
    Blocks are not taken into account here.
    */
    bool can_join_fragments(Fragment* f1, Fragment* f2) const;

    /** Return if two blocks can be joined.
    The blocks must be joinable (Block::can_join) and non empty.
    Both fragments and blocks are taken into account here.
    */
    bool can_join_blocks(Block* b1, Block* b2) const;

    /** Try to join, return empty pointer if failed */
    Block* try_join(Block* one, Block* another) const;

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

    const char* name_impl() const;

private:
    int max_dist_;
    float ratio_to_fragment_;
    float gap_ratio_;
};

}

#endif

