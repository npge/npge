/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FILE_COPY_HPP_
#define BR_FILE_COPY_HPP_

#include "Processor.hpp"

namespace npge {

/** Copy file */
class FileCopy : public Processor {
public:
    /** Constructor */
    FileCopy();

protected:
    void run_impl() const;

    const char* name_impl() const;
};

}

#endif

