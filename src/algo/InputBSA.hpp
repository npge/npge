/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_INPUT_BLOCK_SET_ALIGNMENT_HPP_
#define NPGE_INPUT_BLOCK_SET_ALIGNMENT_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"

namespace npge {

/** Input blockset alignment */
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

