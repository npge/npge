/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FILE_COPY_HPP_
#define BR_FILE_COPY_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Copy file */
class FileCopy : public Processor {
public:
    /** Constructor */
    FileCopy();

protected:
    bool run_impl() const;
};

}

#endif

