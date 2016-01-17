/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PRINT_PARTITION_HPP_
#define NPGE_PRINT_PARTITION_HPP_

#include "AbstractOutput.hpp"

namespace npge {

/** Print overlaps of fragments from target blockset with other as table */
class PrintPartition : public AbstractOutput {
public:
    /** Constructor */
    PrintPartition();

    /** Destructor */
    ~PrintPartition();

protected:
    const char* name_impl() const;

    void prepare() const;

    void print_header(std::ostream& o) const;

    void print_block(std::ostream& o, Block* block) const;

private:
    class Impl;

    Impl* impl_;
};

}

#endif

