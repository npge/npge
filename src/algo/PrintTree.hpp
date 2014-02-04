/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PRINT_TREE_HPP_
#define BR_PRINT_TREE_HPP_

#include "AbstractOutput.hpp"
#include "global.hpp"

namespace bloomrepeats {

/** Print tree.

FragmentDistance is used to calculate distance between fragments.
*/
class PrintTree : public AbstractOutput {
public:
    /** Constructor */
    PrintTree();

    /** Make tree.
    \param block Block.
    \param method Method of tree construction. "upgma" or "nj".
    */
    Tree* make_tree(const Block* block,
            const std::string& method) const;

private:
    FragmentDistance* distance_;

    /** Print table block - tree */
    void print_block(std::ostream& o, Block* block) const;
};

}

#endif

