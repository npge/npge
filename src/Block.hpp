/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_HPP_
#define BR_BLOCK_HPP_

#include <iosfwd>
#include <new>
#include <vector>
#include <string>

#include "global.hpp"

namespace bloomrepeats {

/** Container for fragments.
A block is aimed to keep related fragments together.
*/
class Block {
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

    /** Return a pointer to new instance of Block */
    static BlockPtr create_new();

    /** Constructor.
    Set auto-increment name.
    */
    Block();

    /** Destructor.
    Clear the block.
    */
    ~Block();

    /** Allocate storage */
    void* operator new(size_t x);

    /** Deallocate storage */
    void operator delete(void* ptr);

    /** Return a copy of this block.
    Fragments are copied, sequences are not copied.
    Connections between the fragments
    (\ref Fragment::prev() "prev", \ref Fragment::next() "next")
    are not copied.
    \see BlockSet::clone()
    */
    BlockPtr clone() const;

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
    static int match(BlockPtr one, BlockPtr another);

    /** Filter out and disconnect short and invalid fragments */
    void filter(int min_fragment_length = 100);

    /** Return whether blocks can be joined.
    \param one Block.
    \param another Block.

    Return
     - 1, if fragments of 'one' should preceed fragments from 'another';
     - -1, if fragments of 'another' should preceed fragments from 'one';
     - 0, if blocks can't be joined.

    Blocks can be joined, if they match and all the fragments from
    the first block has an unique neighbor with the same ori
    from the second block.

    \see JoinApprover
    */
    static int can_join(BlockPtr one, BlockPtr another);

    /** Return joined blocks, if these two blocks can be joined.
    Fragments are also \ref Fragment::join "joined".
    \see can_join() for \p logical_ori description.
    */
    static BlockPtr join(BlockPtr one, BlockPtr another, int logical_ori);

    /** Try to join, return empty pointer if failed.
    \param one Block.
    \param another Block.
    \param join_approver Object confirming join.
        Value 0 means always approving one.
    */
    static BlockPtr try_join(BlockPtr one, BlockPtr another,
                             JoinApprover* join_approver = 0);

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

    /** Rearrange this block's fragments before or after neighbors.
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

    This methods adds to the block new fragments, made from neighbor blocks,
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

    /** Return name of block.
    By default, name is auto-increment integer.
    */
    const std::string& name() const {
        return name_;
    }

    /** Set block name.
    Name MUST match the following regex: "[a-zA-Z0-9]{1,40}"
    */
    void set_name(const std::string& name);

private:
    Impl fragments_;
    std::string name_;

    void expand_end(PairAligner& aligner, int batch, int max_overlap);
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Block& block);

}

#endif

