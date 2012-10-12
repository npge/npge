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

    /** Get block set, from which blocks will be added to block_set() */
    const BlockSetPtr& other() const {
        return other_;
    }

    /** Set block set, from which blocks will be added to block_set() */
    void set_other(const BlockSetPtr& other) {
        other_ = other;
    }

private:
    BlockSetPtr other_;
};

}

#endif

