/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_HPP_
#define BR_BLOCK_HPP_

#include <iosfwd>
#include <vector>
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
    typedef std::vector<FragmentPtr> Impl;

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

    /** Return if ori of fragments of two block correspond.
    0 means "no match", 1 means "match as is", -1 means "A match B.inverse()"
    */
    static int match(const BlockPtr& one, const BlockPtr& another);

    /** Return whether blocks can be merged.
    Return
     - 1, if fragments of 'one' should preceed fragments from 'another';
     - -1, if fragments of 'another' should preceed fragments from 'one';
     - 0, if blocks can't be merged.

    Blocks can be merged, if they match and all the fragments from
    the first block has an unique neighbour with the same ori
    from the second block.
    */
    static int can_merge(BlockPtr one, BlockPtr another);

    /** Return merged blocks, if these two blocks can be merged.
    Fragments are also \ref Fragment::merge "merged".
    \see can_merge() for \p logical_ori description.
    */
    static BlockPtr merge(BlockPtr one, BlockPtr another, int logical_ori);

    /** Try to merge, return empty pointer if failed */
    static BlockPtr try_merge(BlockPtr one, BlockPtr another);

    /** Inverse all fragments of this block */
    void inverse();

    /** Max valid shift of the block's fragments.
    Return max value, that can be passed to Fragment::shift_end()
    of each fragment, keeping the fragment Fragment::valid().
    May be negative, if a fragment is already invalid.
    */
    int max_shift_end() const;

    /** Expand block.
    \param aligner Pointer to PairAligner. If aligner = 0,
        then thread specific static one is used.
    \param batch Length of piece, passed to PairAligner at a time.
    \param ori Direction of expansion. 0 means both.
     - One fragment is selected as main.
     - On each iteration, other fragments are aligned to main one.
     - If at least one fragment was aligned on less then 0.5 of batch,
       expansion is stopped.
    */
    void expand(PairAligner* aligner = 0, int batch = 100, int ori = 0);

private:
    Impl fragments_;

    void expand_end(PairAligner& aligner, int batch);

    Block(); // nonconstructible

    friend boost::shared_ptr<Block> boost::make_shared<Block>();
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Block& block);

}

#endif

