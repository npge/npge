/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_INPUT_BLOCK_SET_ALIGNMENT_HPP_
#define BR_INPUT_BLOCK_SET_ALIGNMENT_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"

namespace bloomrepeats {

/** Input block set alignment */
class InputBlockSetAlignment : public Processor {
public:
    /** Constructor */
    InputBlockSetAlignment();

protected:
    bool run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
};

}

#endif

