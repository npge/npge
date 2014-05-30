/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_HASH_HPP_
#define BR_HASH_HPP_

#include "Processor.hpp"
#include "FileWriter.hpp"

namespace npge {

/** Writes hash of blockset */
class Hash : public Processor {
public:
    /** Constructor */
    Hash();

protected:
    void run_impl() const;
    const char* name_impl() const;

private:
    FileWriter file_writer_;
};

}

#endif

