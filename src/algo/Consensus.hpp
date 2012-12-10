/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONSENSUS_HPP_
#define BR_CONSENSUS_HPP_

#include "Processor.hpp"
#include "FileWriter.hpp"

namespace bloomrepeats {

/** Write consensuses of all blocks to a file.
\see ExternalAligner
*/
class Consensus : public Processor, public FileWriter {
protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

