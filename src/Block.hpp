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

    /** Destructor.
    Clear the block.
    */
    virtual ~Block();

    /** Add fragment.
    \attention Two equal fragments must not be inserted!
        For debug build, this is checked with BOOST_ASSERT
    */
    void insert(FragmentPtr fragment);

    /** Remove fragment.
    The fragment is deleted.
    */
    void erase(FragmentPtr fragment);

    /** Return the number of fragments in block */
    size_t size() const;

    /** Return if the block has no fragments */
    bool empty() const;

    /** Return if the block has the fragment */
    bool has(FragmentPtr fragment) const;

    /** Remove all fragments from the block.
    Removed fragments are deleted.
    */
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

    /** Return proportion of columns, composed of size() equal letters.
    \note This method compares fragments as is, without attempt of alignment.
    */
    float identity() const;

    /** Return if ori of fragments of two block correspond.
    0 means "no match", 1 means "match as is", -1 means "A match B.inverse()"
    */
    static int match(const BlockPtr& one, const BlockPtr& another);

    /** Filter out and disconnect short and invalid fragments */
    void filter(int min_fragment_length = 100);

    /** Return whether blocks can be joined.
    \param one Block.
    \param another Block.
    \param max_gap Max \ref Fragment::dist_to "distance between" fragments.

    Return
     - 1, if fragments of 'one' should preceed fragments from 'another';
     - -1, if fragments of 'another' should preceed fragments from 'one';
     - 0, if blocks can't be joined.

    Blocks can be joined, if they match and all the fragments from
    the first block has an unique neighbour with the same ori
    from the second block.
    */
    static int can_join(BlockPtr one, BlockPtr another, size_t max_gap = -1);

    /** Return joined blocks, if these two blocks can be joined.
    Fragments are also \ref Fragment::join "joined".
    \see can_join() for \p logical_ori description.
    */
    static BlockPtr join(BlockPtr one, BlockPtr another, int logical_ori);

    /** Try to join, return empty pointer if failed.
    \param one Block.
    \param another Block.
    \param max_gap Max \ref Fragment::dist_to "distance between" fragments.
    */
    static BlockPtr try_join(BlockPtr one, BlockPtr another,
                             size_t max_gap = -1);

    /** Inverse all fragments of this block.
    \see Fragment::inverse()
    */
    void inverse();

    /** Apply a patch to each fragment of this block.
    \see Fragment::patch()
    */
    void patch(const FragmentDiff& diff);

    /** Split this block into two blocks.
    \see Fragment::split()
    Result must not be null pointer, but may be empty.
    */
    BlockPtr split(size_t new_length);

    /** Rearrange this block's fragments before or after neighbours.
    \see Fragment::find_place()
    */
    void find_place();

    /** Max valid shift of the block's fragments.
    \param max_overlap Max number of positions, that are allowed to be added
       to the block after first overlap occured.
       -1 means "overlaps of any length are allowed".
       Fragments must be \ref BlockSet::connect_fragments "connected"
       for this to work correctly.

    Return max value, that can be passed to Fragment::shift_end()
    of each fragment, keeping the fragment Fragment::valid().
    May be negative, if a fragment is already invalid.
    */
    int max_shift_end(int max_overlap = 0) const;

    /** Expand block.
    \param aligner Pointer to PairAligner. If aligner = 0,
        then thread specific static one is used.
    \param batch Length of piece, passed to PairAligner at a time.
    \param ori Direction of expansion. 0 means both.
    \param max_overlap Max number of positions, that are allowed to be added
       to the block after first overlap occured.
       -1 means "overlaps of any length are allowed".
       Fragments must be \ref BlockSet::connect_fragments "connected"
       for this to work correctly.

    Steps:
     - One fragment is selected as main.
     - On each iteration, other fragments are aligned to main one.
     - If at least one fragment was aligned on less then 0.5 of batch,
       expansion is stopped.
    */
    void expand(PairAligner* aligner = 0, int batch = 100, int ori = 0,
                int max_overlap = 0);

    /** Return number of the fragment's positions, occupied by the block */
    size_t common_positions(const Fragment& fragment);

    /** Expand block by fragments.
    \param aligner Pointer to PairAligner. If aligner = 0,
        then thread specific static one is used.
    \param batch Length of piece, passed to PairAligner at a time.

    This methods adds to the block new fragments, made from neighbour blocks,
    if they are \ref PairAligner::aligned() "aligned" with
    some fragment from this block.

    Return \p true is something was added.

    \warning
       Fragments must be \ref BlockSet::connect_fragments "connected"
       for this to work correctly.
    */
    bool expand_by_fragments(PairAligner* aligner = 0, int batch = 100);

    /** Move contents of other to this.
    Other is cleared.
    Duplicates are removed (\p other is \ref inverse "inversed" if needed).
    */
    void merge(BlockPtr other);

private:
    Impl fragments_;

    void expand_end(PairAligner& aligner, int batch, int max_overlap);

    Block(); // nonconstructible

    friend BlockPtr boost::make_shared<Block>();
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Block& block);

}

#endif

