/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_SET_HPP_
#define BR_BLOCK_SET_HPP_

#include <iosfwd>
#include <set>

#include "global.hpp"

namespace bloomrepeats {

/** Container of blocks.
*/
class BlockSet {
public:
    /** Type of implementation container.
    Do not rely on ths type!
    */
    typedef std::set<BlockPtr> Impl;

    /** Iterator */
    typedef Impl::iterator iterator;

    /** Constant iterator */
    typedef Impl::const_iterator const_iterator;

    /** Add block.
    The same block can't be added twice.
    */
    void insert(BlockPtr block);

    /** Remove fragment */
    void erase(BlockPtr block);

    /** Return the number of blocks */
    size_t size() const;

    /** Return if there is no blocks */
    bool empty() const;

    /** Return if has the block */
    bool has(BlockPtr block) const;

    /** Remove all blocks */
    void clear();

    /** Get some block if any or an empty pointer */
    BlockPtr front() const;

    /** Return iterator to beginning */
    iterator begin();

    /** Return constant iterator to beginning */
    const_iterator begin() const;

    /** Return iterator to end */
    iterator end();

    /** Return constant iterator to end */
    const_iterator end() const;

    /** Connect all the fragments (prev-next) */
    void connect_fragments();

    /** Filter out fragments and blocks */
    void filter(int min_fragment_length = 100, int min_block_size = 2);

    /** Merge neighbour blocks */
    void merge();

    /** Expand all blocks (starting from blocks of large number of fragments).
    \see Block::expand()
    */
    void expand_blocks(PairAligner* aligner = 0, int batch = 100, int ori = 0,
                       bool overlap = false);

private:
    Impl blocks_;
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const BlockSet& block_set);

}

#endif

