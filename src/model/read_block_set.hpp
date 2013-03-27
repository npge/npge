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

#include "global.hpp"
#include "FastaReader.hpp"

namespace bloomrepeats {

/** Read fasta representation of blocks to block set */
class BlockSetFastaReader : public FastaReader {
public:
    /** Constructor.
    \param block_set BlockSet to read to.
    \param input Input stream.
    \param keep_alignment If the alignment is read too.
    \param type Storage type of alignment rows.
    \param seq_type If not NO_SEQUENCE, sequences are read too, of this type.
        Blocks must cover whole sequences.
        Sequences can be passed as is, before fragments, then effect is
        like from AddSequences.
    */
    BlockSetFastaReader(BlockSet& block_set, std::istream& input,
                        bool keep_alignment, RowType type,
                        SequenceType seq_type = NO_SEQUENCE);

protected:
    void new_sequence(const std::string& name, const std::string& description);

    void grow_sequence(const std::string& data);

private:
    BlockSet& block_set_;
    std::map<std::string, Block*> name2block_;
    bool keep_alignment_;
    RowType row_type_;
    SequenceType seq_type_;
    AlignmentRow* row_;
    Fragment* fragment_; // 0 if read whole sequence
    Sequence* sequence_; // not 0, if read whole seq or seq from fragment
    size_t used_np_;
};

}

#endif

