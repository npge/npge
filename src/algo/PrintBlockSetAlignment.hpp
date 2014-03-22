/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PRINT_BLOCK_SET_ALIGNMENT_HPP_
#define BR_PRINT_BLOCK_SET_ALIGNMENT_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "FileWriter.hpp"

namespace bloomrepeats {

/** Print block set alignment */
class PrintBlockSetAlignment : public Processor {
public:
    /** Constructor */
    PrintBlockSetAlignment();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileWriter file_writer_;
};

}

#endif

