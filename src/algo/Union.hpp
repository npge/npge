/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_UNION_HPP_
#define NPGE_UNION_HPP_

#include "Processor.hpp"

namespace npge {

/** Add clones of blocks from another blockset to this blockset.
Deep copy is performed. Alignment rows are copied.
*/
class Union : public Processor {
public:
    /** Constructor
    \param source BlockSet, from which blocks will be added to block_set().
    */
    Union(const BlockSetPtr& source = BlockSetPtr());

    /** Return a copy of this fragment.
    Alignment row is copied.
    */
    static Fragment* clone_fragment(Fragment* source);

    /** Return a copy of this block.
    Fragments are copied, sequences are not copied.
    */
    static Block* clone_block(Block* source);

    /** Return a copy of this blockset.
    Fragments and blocks are copied, sequences are not copied,
    sequence list is copied.
    */
    static BlockSetPtr clone_block_set(BlockSetPtr block_set);

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

