/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ALL_OPTIONS_HPP_
#define NPGE_ALL_OPTIONS_HPP_

#include "Processor.hpp"
#include "FileWriter.hpp"

namespace npge {

/** Print HTML help about all global configuration options */
class AllOptions : public Processor {
public:
    /** Constructor */
    AllOptions();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileWriter out_;
};

}

#endif

