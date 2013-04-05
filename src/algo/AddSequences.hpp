/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ADD_SEQUENCES_HPP_
#define BR_ADD_SEQUENCES_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"
#include "SeqStorage.hpp"

namespace bloomrepeats {

/** Add input sequences to the block set.
\deprecated Use AddBlocks
*/
class AddSequences : public Processor, public FileReader, public SeqStorage {
public:
    /** Constructor */
    AddSequences(SequenceType seq_type = COMPACT_SEQUENCE);

protected:
    /** Add options to options description */
    void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

