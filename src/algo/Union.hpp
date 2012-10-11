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

/** Add clones of blocks from another block set to this block set */
class Union : public Processor {
public:
    /** Constructor
    \param source BlockSet, from which blocks will be added to block_set().
    */
    Union(const BlockSetPtr& source);

    /** Get block set, from which blocks will be added to block_set() */
    const BlockSetPtr& source() const {
        return source_;
    }

    /** Set block set, from which blocks will be added to block_set() */
    void set_source(const BlockSetPtr& source) {
        source_ = source;
    }

protected:
    /** Apply the action */
    bool run_impl() const;

private:
    BlockSetPtr source_;
};

}

#endif

