/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PRINT_OVERLAPS_HPP_
#define NPGE_PRINT_OVERLAPS_HPP_

#include <vector>

#include "AbstractOutput.hpp"
#include "FragmentCollection.hpp"

namespace npge {

/** Print ASCII diagram with all fragments overlapping with a block.

It is recommended to use this processor if blocks have alignment.
*/
class PrintOverlaps : public AbstractOutput {
public:
    /** Constructor */
    PrintOverlaps();

    /** Print ASCII diagram with all fragments overlapping with a block */
    void print_block(std::ostream& o, Block* block) const;

protected:
    const char* name_impl() const;
    void prepare() const;
    void finish_work_impl() const;

private:
    mutable VectorFc s2f_; // FIXME
};

}

#endif

