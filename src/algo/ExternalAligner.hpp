/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_EXTERNAL_ALIGNER_HPP_
#define BR_EXTERNAL_ALIGNER_HPP_

#include "BlocksJobs.hpp"

namespace bloomrepeats {

/** Align blocks with external alignment tool.
Skips block, if block's fragment has row.
*/
class ExternalAligner : public BlocksJobs {
public:
    /** Constructor
    \param cmd Command template. Use %1% as input of aligner, %2% as output.
    */
    ExternalAligner(const std::string& cmd =
                        "mafft --quiet --retree 1 --maxiterate 1 %1% > %2%");

    /** Apply external aligner to a blick */
    void align_block(Block* block) const;

protected:
    void change_blocks_impl(std::vector<Block*>& blocks) const;

    void process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;
};

}

#endif

