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

#include "global.hpp"

namespace npge {

/** Read fasta representation of blocks to block set */
class BlockSetFastaReader {
public:
    /** Constructor.
    \param block_set BlockSet to read to ("target").
    \param input Input stream.
    \param type Storage type of alignment rows.
    \param seq_type Type of sequences created by this processor.
    */
    BlockSetFastaReader(BlockSet& block_set, std::istream& input,
                        RowType type, SequenceType seq_type);

    /** Destructor */
    virtual ~BlockSetFastaReader();

    /** Associate name with block set */
    void set_block_set(const std::string& name,
                       BlockSet* block_set);

    /** Return blockset by name or 0 */
    BlockSet* get_block_set(const std::string& name) const;

    /** Return if unknown block set name is silently skipped */
    bool unknown_bs_allowed() const;

    /** Set if unknown block set name is silently skipped.
    Otherwise Exception is thrown.
    Defaults to true.
    */
    void set_unknown_bs_allowed(bool unknown_bs_allowed);

    /** Return number of workers */
    int workers() const;

    /** Set number of workers */
    void set_workers(int workers);

    /** Run the reader */
    void run();

private:
    class Impl;
    Impl* impl_;
};

}

#endif

