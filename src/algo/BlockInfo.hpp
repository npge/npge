/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_BLOCK_INFO_HPP_
#define NPGE_BLOCK_INFO_HPP_

#include "AbstractOutput.hpp"

namespace npge {

/** Print information about blocks.

It is recommended to use this processor if blocks have alignment.
*/
class BlockInfo : public AbstractOutput {
public:
    /** Constructor */
    BlockInfo(const std::string& prefix = "info-");

    /** Print information about block */
    void print_block(std::ostream& o, Block* block) const;

    void print_header(std::ostream& o) const;

protected:
    const char* name_impl() const;
};

}

#endif

