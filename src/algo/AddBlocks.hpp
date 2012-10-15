/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ADD_BLOCKS_HPP_
#define BR_ADD_BLOCKS_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "FileReader.hpp"

namespace bloomrepeats {

/** Add input blocks to the block set.
\note This processor depends on AddSequences.

Wrapper for stream >> block_set or stream >> alignment.
*/
class AddBlocks : public Processor, public FileReader {
public:
    /** Default constructor.
    Read block set only.
    */
    AddBlocks();

    /** Default constructor.
    Read block set and alignment.
    If alignment.block_set() != block_set,
    then alignment.set_block_set(block_set()) is called.
    */
    AddBlocks(const AlignmentPtr& alignment);

    /** Get alignment associated with the block set */
    const AlignmentPtr& alignment() const {
        return alignment_;
    }

    /** Set alignment associated with the block set */
    void set_alignment(const AlignmentPtr& alignment) {
        alignment_ = alignment;
    }

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

    const char* name_impl() const;

private:
    AlignmentPtr alignment_;
};

}

#endif

