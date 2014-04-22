/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_EXTERNAL_ALIGNER_HPP_
#define BR_EXTERNAL_ALIGNER_HPP_

#include <vector>

#include "BlocksJobs.hpp"
#include "global.hpp"

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

    /** Apply external aligner to a block */
    void align_block(Block* block) const;

    /** Apply sequences */
    void align_seqs(Strings& seqs) const;

    /** Return if alignment is needed and build it in obvious cases */
    bool alignment_needed(Block* block) const;

    /** Apply external aligner to file */
    void align_file(const std::string& input,
                    const std::string& output) const;

    /** Return list of alignment rows from fasta file */
    void read_alignment(Strings& rows,
                        const std::string& file) const;

protected:
    void change_blocks_impl(std::vector<Block*>& blocks) const;

    void process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;
};

}

#endif

