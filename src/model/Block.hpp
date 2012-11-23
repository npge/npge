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

/** Column stat of alignment */
struct AlignmentStat {
    /** Default constructor */
    AlignmentStat();

    /** Non-empty ident columns without gaps */
    int ident_nogap;

    /** Non-empty ident columns with gaps */
    int ident_gap;

    /** Non-empty non-ident columns without gaps */
    int noident_nogap;

    /** Non-empty non-ident columns with gaps */
    int noident_gap;

    /** Empty columns (consist only of gaps) */
    int pure_gap;

    /** All columns.
    Must be equal to sum of above variables.
    */
    int total;
};

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

    /** Return length of alignment.
    If a fragment doesn't have alignment row attached,
    then length of the fragment is taken.
    */
    size_t alignment_length() const;

    /** Make alignment stat of alignment */
    void make_stat(AlignmentStat& stat) const;

    /** Return proportion of columns, composed of size() equal letters.
    If a fragment doesn't have alignment row attached,
    then it is taken as is.
    */
    float identity() const;

    /** Return consensus letter for given position.
    For each column, the most frequent letter is written to consensus.
    If frequencies of several letters are equal, them some of them is written.
    For pure gap columns, value of argument 'gap' is written.
    */
    char consensus_char(int pos, char gap = 'a') const;

    /** Write consensus to output stream.
    \see consensus_char
    */
    void consensus(std::ostream& o, char gap = 'a') const;

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

private:
    Impl fragments_;
    std::string name_;
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const Block& block);

}

#endif

