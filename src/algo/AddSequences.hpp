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

namespace bloomrepeats {

/** Add input sequences to the block set.
\deprecated Use AddBlocks
*/
class AddSequences : public Processor {
public:
    /** Constructor */
    AddSequences();

protected:
    /** Apply the action */
    bool run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
};

}

#endif

