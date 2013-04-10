/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_UNION_HPP_
#define BR_UNION_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Add clones of blocks from another block set to this block set.
Deep copy is performed. Alignment rows are copied.
*/
class Union : public Processor {
public:
    /** Constructor
    \param source BlockSet, from which blocks will be added to block_set().
    */
    Union(const BlockSetPtr& source = BlockSetPtr());

    /** Return a copy of this block.
    Fragments are copied, sequences are not copied.
    Connections between the fragments
    (\ref Fragment::prev() "prev", \ref Fragment::next() "next")
    are not copied.
    */
    static Block* clone_block(Block* source);

    /** Return a copy of this block set.
    Fragments and blocks are copied, sequences are not copied,
    sequence list is copied.
    Connections between the fragments
    (\ref Fragment::prev() "prev", \ref Fragment::next() "next")
    are rebuild with Connector.
    \todo Preserve fragment connections from source block set.
    */
    static BlockSetPtr clone_block_set(BlockSetPtr block_set);

protected:
    /** Apply the action */
    bool run_impl() const;
};

}

#endif

