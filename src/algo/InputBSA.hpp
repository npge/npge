/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_INPUT_BLOCK_SET_ALIGNMENT_HPP_
#define BR_INPUT_BLOCK_SET_ALIGNMENT_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"

namespace npge {

/** Input block set alignment */
class InputBSA : public Processor {
public:
    /** Constructor */
    InputBSA();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
};

}

#endif

