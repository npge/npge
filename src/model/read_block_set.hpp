/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_READ_BLOCK_SET_HPP_
#define NPGE_READ_BLOCK_SET_HPP_

#include <iosfwd>
#include <string>
#include <vector>
#include <map>
#include <set>

#include "global.hpp"
#include "FastaReader.hpp"

namespace npge {

/** Read fasta representation of blocks to block set */
class BlockSetFastaReader : public FastaReader {
public:
    /** Constructor.
    \param block_set BlockSet to read to ("target").
    \param input Input stream.
    \param type Storage type of alignment rows.
    \param seq_type Type of sequences created by this processor.
    */
    BlockSetFastaReader(BlockSet& block_set, std::istream& input,
                        RowType type, SequenceType seq_type);

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
    typedef std::map<std::string, BlockSet*> Name2BlockSet;
    typedef std::map<std::string, Block*> Name2Block;
    typedef std::map<std::string, SequencePtr> Name2Seq;
    typedef std::map<BlockSet*, Name2Block> Bs2Name2Block;

    std::vector<BlockSet*> block_sets_;
    bool unknown_bs_allowed_;
    Name2BlockSet name2block_set_;
    Bs2Name2Block name2block_;
    Name2Seq name2seq_;
    bool keep_alignment_;
    RowType row_type_;
    SequenceType seq_type_;
    std::vector<AlignmentRow*> rows_;
    Fragment* fragment_; // 0 if read whole sequence
    Sequence* sequence_; // not 0, if read whole seq or seq from fragment
    size_t used_np_;
    std::set<Sequence*> seqs_from_frags_; // sequences read from fragments
};

}

#endif

