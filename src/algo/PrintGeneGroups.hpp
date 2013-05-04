/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PRINT_GENE_GROUPS_HPP_
#define BR_PRINT_GENE_GROUPS_HPP_

#include "AbstractOutput.hpp"

namespace bloomrepeats {

/** Print information about gene groups.
Block sets: 'target' (gene groups, weak blocks), 'pangenome'.
*/
class PrintGeneGroups : public AbstractOutput {
public:
    /** Constructor */
    PrintGeneGroups();

    /** Destructor */
    ~PrintGeneGroups();

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

