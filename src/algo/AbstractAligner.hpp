/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ABSTRACT_ALIGNER_HPP_
#define NPGE_ABSTRACT_ALIGNER_HPP_

#include <vector>

#include "BlocksJobs.hpp"
#include "global.hpp"

namespace npge {

/** Align blocks.
Skips block, if block's fragment has row.
*/
class AbstractAligner : public BlocksJobs {
public:
    /** Constructor */
    AbstractAligner();

    /** Align a block */
    void align_block(Block* block) const;

    /** Apply sequences */
    void align_seqs(Strings& seqs) const;

    /** Return if alignment is needed and build it in obvious cases */
    bool alignment_needed(Block* block) const;

    /** Remove pure gap columns */
    static void remove_pure_gap_columns(Block* block);

    /** Return if the aligner works */
    bool test(bool gaps = false) const;

    /** Return aligner type */
    virtual std::string aligner_type() const = 0;

protected:
    void change_blocks_impl(Blocks& blocks) const;

    void process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;

    /** Align sequences.
    seqs is guaranteed not to be empty.
    Each sequence is guaranteed not to be empty.
    */
    virtual void align_seqs_impl(Strings& seqs) const = 0;
};

}

#endif

