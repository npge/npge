/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_HPP_
#define BR_BLOCK_HPP_

#include <set>
#include <boost/enable_shared_from_this.hpp>

#include "global.hpp"

namespace bloomrepeats {

/** Container for fragments.
A block is aimed to keep related fragments together.
*/
class Block : public boost::enable_shared_from_this<Block> {
public:
    /** Type of implementation container.
    Do not rely on ths type!
    To traverse all fragments, use BOOST_FOREACH (FragmentPtr f, block).
    For other operations use public members of Block.
    */
    typedef std::set<FragmentPtr> Impl;

    /** Iterator */
    typedef Impl::iterator iterator;

    /** Constant iterator */
    typedef Impl::const_iterator const_iterator;

    /** Instead of constructor */
    static BlockPtr create_new();

    /** Add fragment.
    \attention Two equal fragments must not be inserted!
        For debug build, this is checked with BOOST_ASSERT
    */
    void insert(FragmentPtr fragment);

    /** Remove fragment */
    void erase(FragmentPtr fragment);

    /** Return the number of fragments in block */
    size_t size() const;

    /** Return if the block has no fragments */
    bool empty() const;

    /** Return if the block has the fragment */
    bool has(FragmentPtr fragment) const;

    /** Remove all fragments from the block */
    void clear();

    /** Get some fragment if any or an empty pointer */
    FragmentPtr front() const;

    /** Return iterator to beginning */
    iterator begin();

    /** Return constant iterator to beginning */
    const_iterator begin() const;

    /** Return iterator to end */
    iterator end();

    /** Return constant iterator to end */
    const_iterator end() const;

    /** Inverse all fragments of this block */
    void inverse();

    /** Expand block.
    \param batch Length of piece, passed to PairAligner at a time.
    \param max_errors Max number of errors in pair alignment,
        that may happen on a batch.
    \param gap_range Max distance from main diagonal of considered
        states of pair alignment. The more gap_range, the more time.
    \param ori Direction of expansion. 0 means both.
     - One fragment is selected as main.
     - On each iteration, other fragments are aligmned to main one.
     - If at least one fragment was aligned on less then 0.5 of batch,
       expansion is interrupted (at previous step).
    */
    void expand(int batch = 100, int max_errors = 5, int gap_range = 5,
                int ori = 0);

private:
    std::set<FragmentPtr> fragments_;

    void expand_end(int batch, int max_errors, int gap_range);

    Block(); // nonconstructible

    friend boost::shared_ptr<Block> boost::make_shared<Block>();
};

}

#endif

