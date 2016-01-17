/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OUTPUT_PIPE_HPP_
#define NPGE_OUTPUT_PIPE_HPP_

#include "Pipe.hpp"

namespace npge {

/** Grace output */
class Write : public Pipe {
public:
    /** Constructor */
    Write(const std::string& prefix = "out-");

protected:
    const char* name_impl() const;
};

}

#endif

