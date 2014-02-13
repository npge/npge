/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PRINT_TREE_HPP_
#define BR_PRINT_TREE_HPP_

#include "AbstractOutput.hpp"
#include "tree.hpp"
#include "global.hpp"

namespace bloomrepeats {

class FragmentLeaf : public LeafNode {
public:
    FragmentLeaf(const Fragment* f, const FragmentDistance* distance);

    double distance_to_impl(const LeafNode* leaf) const;

    std::string name_impl() const;

    TreeNode* clone_impl() const;

    const Fragment* fragment() const {
        return f_;
    }

private:
    const Fragment* f_;
    const FragmentDistance* distance_;
};

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
    TreeNode* make_tree(const Block* block,
                        const std::string& method) const;

    /** Make tree.
    \param block Block.
    */
    TreeNode* make_tree(const Block* block) const;

private:
    FragmentDistance* distance_;

    /** Print table block - tree */
    void print_block(std::ostream& o, Block* block) const;

    void print_header(std::ostream& o) const;
};

}

#endif

