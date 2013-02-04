/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OTHER_BLOCK_SET_HPP_
#define BR_OTHER_BLOCK_SET_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Add clones of blocks from another block set to this block set */
class OtherBlockSet {
public:
    /** Constructor
    \param other BlockSet.
    */
    OtherBlockSet(const BlockSetPtr& other);

    /** Get other set.
    If OtherBlockSet was initialized with a ProcessorPtr and set_other was not
    called, then return current block_set() of that Processor.
    */
    BlockSetPtr other() const;

    /** Set other block set */
    void set_other(const BlockSetPtr& other);

    /** Set other block set */
    void set_processor(const Processor* processor);

    /** Set other block set */
    void set_other_block_set(OtherBlockSet* other_block_set);

private:
    BlockSetPtr other_;
    const Processor* processor_;
    OtherBlockSet* other_block_set_;
};

}

#endif

