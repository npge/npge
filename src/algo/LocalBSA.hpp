/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_LOCAL_BSA_HPP_
#define NPGE_LOCAL_BSA_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace npge {

class FindBSA;

/** Build blockset alignments for global blocks */
class LocalBSA : public Processor {
public:
    /** Constructor */
    LocalBSA();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FindBSA* bsa_;
};

}

#endif

