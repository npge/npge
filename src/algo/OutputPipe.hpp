/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OUTPUT_PIPE_HPP_
#define NPGE_OUTPUT_PIPE_HPP_

#include "Pipe.hpp"

namespace npge {

/** OriByMajority, Connector, UniqueNames, Output */
class OutputPipe : public Pipe {
public:
    /** Constructor */
    OutputPipe(const std::string& prefix = "out-");

protected:
    const char* name_impl() const;
};

}

#endif

