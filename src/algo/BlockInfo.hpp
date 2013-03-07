/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLOCK_INFO_HPP_
#define BR_BLOCK_INFO_HPP_

#include "AbstractOutput.hpp"

namespace bloomrepeats {

/** Print information about blocks.

It is recommended to use this processor if blocks have alignment.
*/
class BlockInfo : public AbstractOutput {
public:
    /** Print information about block */
    void print_block(std::ostream& o, Block* block) const;

protected:
    const char* name_impl() const;
};

}

#endif

