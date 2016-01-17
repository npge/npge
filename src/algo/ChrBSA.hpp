/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_CHR_BLOCK_SET_ALIGNMENT_HPP_
#define NPGE_CHR_BLOCK_SET_ALIGNMENT_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace npge {

class FindBSA;

/** Build blockset alignments for all chromosomes */
class ChrBSA : public Processor {
public:
    /** Constructor */
    ChrBSA();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FindBSA* bsa_;
};

}

#endif

