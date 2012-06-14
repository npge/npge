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

    /** Filter out fragments and blocks.
    \see Block::filter()
    */
    void filter(int min_fragment_length = 100, int min_block_size = 2);

    /** Merge neighbour blocks.
    \param max_gap Max \ref Fragment::dist_to "distance between" fragments.
    */
    void join(size_t max_gap = -1);

    /** Expand all blocks (starting from blocks of large number of fragments).
    \see Block::expand()
    */
    void expand_blocks(PairAligner* aligner = 0, int batch = 100, int ori = 0,
                       bool overlap = false);

    /** Return if there are blocks which have intersecting fragments.
    \warning Fragments must be \ref BlockSet::connect_fragments "connected"
    */
    bool intersections() const;

    /** Resolve intersecting fragments.
    If some blocks from the block set have intersecting fragments
    these two blocks are replaced with one higher (and narrower)
    block and several remainder blocks.

    Anyway, applying this method guarantees that no blocks of the block set
    have intersecting fragments.

    Since resolve_intersections() can split blocks,
    applying \ref join "join(0)" is recommended.

    \verbatim
    Input:
        Block 1:
            seq1: ---xxxx----
            seq2: ---xxxx----
        Block 2
            seq2: -----xxxx--
            seq3: -----xxxx--
            seq4: -----xxxx--

    Output of resolve_intersections:
        Block 1:
            seq1: ---xx------
            seq2: ---xx------
        Block 2
            seq2: -------xx--
            seq3: -------xx--
            seq4: -------xx--
        Block 3:
            seq1: -----xx----
            seq2: -----xx----
            seq3: -----xx----
            seq4: -----xx----
    \endverbatim
    \warning Fragments must be \ref BlockSet::connect_fragments "connected"
       for this to work correctly.
    */
    void resolve_intersections();

    /** Expand all blocks by fragments.
    Return \p true is something was added.
    \see Block::expand_by_fragments().
    */
    bool expand_blocks_by_fragments(PairAligner* aligner = 0, int batch = 100);

private:
    Impl blocks_;

    BlockPtr treat_two(const FragmentPtr& x, const FragmentPtr& y,
                       int min_intersection);
};

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const BlockSet& block_set);

}

#endif

