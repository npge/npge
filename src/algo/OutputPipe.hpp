/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OUTPUT_PIPE_HPP_
#define BR_OUTPUT_PIPE_HPP_

#include "Pipe.hpp"

namespace bloomrepeats {

/** OriByMajority, Connector, UniqueNames, Output */
class OutputPipe : public Pipe {
public:
    /** Constructor */
    OutputPipe(const std::string& prefix = "out-");
};

}

#endif

