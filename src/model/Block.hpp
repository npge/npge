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
    To traverse all fragments, use BOOST_FOREACH (Fragment* f, block).
    For other operations use public members of Block.
    */
    typedef std::vector<Fragment*> Impl;

    /** Iterator */
    typedef Impl::iterator iterator;

    /** Constant iterator */
    typedef Impl::const_iterator const_iterator;

    /** Constructor.
    Set random name.
    */
    Block();

    /** Constructor */
    Block(const std::string& name);

    /** Destructor.
    Clear the block.
    */
    ~Block();

    /** Allocate storage */
    void* operator new(size_t x);

    /** Deallocate storage */
    void operator delete(void* ptr);

    /** Add fragment.
    \attention Two equal fragments must not be inserted!
        For debug build, this is checked with BOOST_ASSERT
    */
    void insert(Fragment* fragment);

    /** Remove fragment.
    The fragment is deleted.
    */
    void erase(Fragment* fragment);

    /** Return the number of fragments in block */
    size_t size() const;

    /** Return if the block has no fragments */
    bool empty() const;

    /** Return if the block has the fragment */
    bool has(Fragment* fragment) const;

    /** Remove all fragments from the block.
    Removed fragments are deleted.
    */
    void clear();

    /** Exchange values of two objects */
    void swap(Block& other);

    /** Get some fragment if any or an empty pointer */
    Fragment* front() const;

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
    static int match(Block* one, Block* another);

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
    Block* split(size_t new_length);

    /** Rearrange this block's fragments before or after neighbors.
    \see Fragment::find_place()
    */
    void find_place();

    /** Max valid shift of the block's fragments.
    \param max_overlap Max number of positions, that are allowed to be added
       to the block after first overlap occured.
       -1 means "overlaps of any length are allowed".
       Fragments must be \ref Connector "connected"
       for this to work correctly.

    Return max value, that can be passed to Fragment::shift_end()
    of each fragment, keeping the fragment Fragment::valid().
    May be negative, if a fragment is already invalid.
    */
    int max_shift_end(int max_overlap = 0) const;

    /** Return number of the fragment's positions, occupied by the block */
    size_t common_positions(const Fragment& fragment);

    /** Move contents of other to this.
    Other is cleared.
    Duplicates are removed (\p other is \ref inverse "inversed" if needed).
    */
    void merge(Block* other);

    /** Return name of block.
    By default, name is "00000000".
    */
    const std::string& name() const {
        return name_;
    }

    /** Set block name.
    Name MUST match the following regex: "[a-zA-Z0-9]{1,40}"
    Name MUST be unique (it is not checked by this method).
    */
    void set_name(const std::string& name);

    /** Set random name */
    void set_random_name();

    /** Set name based on fragments.
    This name will be same for same sets on child fragments.
    */
    void set_name_from_fragments();

    /** Get alignment of this block.
    If alignment was not set, this method returns 0.
    */
    Alignment* alignment() const {
        return alignment_;
    }

    /** Set alignment of this block.
    By default, alignment is not set.

    Ownership if transferred.
    Previous alignment of this block, if present, is deleted.
    If the alignment passed is already owned by a block,
    it is disconnected from it.

    alignment may be 0.
    */
    void set_alignment(Alignment* alignment);

private:
    Impl fragments_;
    std::string name_;
    Alignment* alignment_;

    friend class Alignment;
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Block& block);

}

#endif

