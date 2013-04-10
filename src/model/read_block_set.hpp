/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_READ_BLOCK_SET_HPP_
#define BR_READ_BLOCK_SET_HPP_

#include <iosfwd>
#include <string>
#include <map>
#include <set>

#include "global.hpp"
#include "FastaReader.hpp"

namespace bloomrepeats {

/** Read fasta representation of blocks to block set */
class BlockSetFastaReader : public FastaReader {
public:
    /** Constructor.
    \param block_set BlockSet to read to ("target").
    \param input Input stream.
    \param keep_alignment If the alignment is read too.
    \param type Storage type of alignment rows.
    \param seq_type Type of sequences created by this processor.
    */
    BlockSetFastaReader(BlockSet& block_set, std::istream& input,
                        bool keep_alignment, RowType type,
                        SequenceType seq_type);

    /** Associate name with block set */
    void set_block_set(const std::string& name, BlockSet* block_set);

    /** Return if unknown block set name is silently skipped */
    bool unknown_bs_allowed() const {
        return unknown_bs_allowed_;
    }

    /** Set if unknown block set name is silently skipped.
    Otherwise Exception is thrown.
    Defaults to true.
    */
    void set_unknown_bs_allowed(bool unknown_bs_allowed) {
        unknown_bs_allowed_ = unknown_bs_allowed;
    }

protected:
    void new_sequence(const std::string& name, const std::string& description);

    void grow_sequence(const std::string& data);

private:
    BlockSet* block_set_;
    bool unknown_bs_allowed_;
    std::map<std::string, BlockSet*> name2block_set_;
    std::map<std::string, Block*> name2block_;
    bool keep_alignment_;
    RowType row_type_;
    SequenceType seq_type_;
    AlignmentRow* row_;
    Fragment* fragment_; // 0 if read whole sequence
    Sequence* sequence_; // not 0, if read whole seq or seq from fragment
    size_t used_np_;
    std::set<Sequence*> seqs_from_frags_; // sequences read from fragments
};

}

#endif

