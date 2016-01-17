/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GLOBAL_TREE_HPP_
#define NPGE_GLOBAL_TREE_HPP_

#include "global.hpp"
#include "Processor.hpp"
#include "FileWriter.hpp"

namespace npge {

/** Print global tree.

Stem blocks are selected.
Distances between fragments are summed,
then NJ-tree is built.
*/
class GlobalTree : public Processor {
public:
    /** Constructor */
    GlobalTree();

protected:
    void run_impl() const;
    const char* name_impl() const;

private:
    FileWriter tre_;
};

}

#endif

