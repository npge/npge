/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FRAGMENTS_EXTENDER_HPP_
#define NPGE_FRAGMENTS_EXTENDER_HPP_

#include "BlocksJobs.hpp"

namespace npge {

class MetaAligner;

/** Move block's boundaries and align only new parts.
Blocks without alignment and blocks of <= 2 fragments are not changed.
*/
class FragmentsExtender : public BlocksJobs {
public:
    /** Constructor */
    FragmentsExtender();

    /** Extend one block */
    void extend(Block* block) const;

protected:
    void process_block_impl(Block* block, ThreadData*) const;
    const char* name_impl() const;

private:
    MetaAligner* aligner_;
};

}

#endif

