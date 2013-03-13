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

Wrapper for stream >> block_set or stream >> alignment_row.
*/
class AddBlocks : public Processor, public FileReader {
public:
    /** Default constructor.
    Read block set only.
    */
    AddBlocks(bool keep_alignment = false);

    /** Get if alignments is extracted too (not only blocks) */
    bool keep_alignment() const {
        return keep_alignment_;
    }

    /** Set if alignments is extracted too (not only blocks) */
    void set_keep_alignment(bool keep_alignment) {
        keep_alignment_ = keep_alignment;
    }

    /** Get alignment row type.
     - "map": MapAlignmentRow
     - "compact": CompactAlignmentRow

    Defaults to 'compact'.
    */
    const std::string& row_type() const {
        return row_type_;
    }

    /** Set alignment row type */
    void set_row_type(const std::string& row_type) {
        row_type_ = row_type;
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
    bool keep_alignment_;
    std::string row_type_;
};

}

#endif

