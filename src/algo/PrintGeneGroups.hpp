/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PRINT_GENE_GROUPS_HPP_
#define NPGE_PRINT_GENE_GROUPS_HPP_

#include "AbstractOutput.hpp"

namespace npge {

/** Print information about gene groups.
Blocksets: 'target' (gene groups, weak blocks), 'pangenome'.
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

