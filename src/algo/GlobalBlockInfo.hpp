/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GLOBAL_BLOCK_INFO_HPP_
#define NPGE_GLOBAL_BLOCK_INFO_HPP_

#include <map>
#include <set>

#include "AbstractOutput.hpp"

namespace npge {

/** Print information about global blocks */
class GlobalBlockInfo : public AbstractOutput {
public:
    /** Constructor */
    GlobalBlockInfo(const std::string& prefix = "ginfo-");

    /** Print information about block */
    void print_block(std::ostream& o, Block* block) const;

    void print_header(std::ostream& o) const;

protected:
    const char* name_impl() const;

    void prepare() const;

    typedef std::set<Block*> BlocksSet;
    typedef std::map<Block*, BlocksSet> G2N;
    mutable G2N global2normal_;
};

}

#endif
