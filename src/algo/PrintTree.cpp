/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/cast.hpp>
#include <boost/foreach.hpp>

#include "PrintTree.hpp"
#include "FragmentDistance.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "Exception.hpp"
#include "tree.hpp"

namespace bloomrepeats {

PrintTree::PrintTree() {
    distance_ = new FragmentDistance;
    distance_->set_parent(this);
    set_opt_prefix("tree-");
    add_opt("method",
            "Method of tree construction (upgma/nj)",
            std::string("nj"));
}

class FragmentLeaf : public LeafNode {
public:
    FragmentLeaf(const Fragment* f, const FragmentDistance* distance):
        f_(f), distance_(distance)
    { }

    float distance_to_impl(const LeafNode* leaf) const {
        const FragmentLeaf* fl;
        fl = boost::polymorphic_downcast<const FragmentLeaf*>(leaf);
        return distance_->fragment_distance(f_, fl->f_).ratio();
    }

    std::string name_impl() const {
        return f_->id();
    }

    AbstractTreeNode* clone_impl() const {
        return new FragmentLeaf(f_, distance_);
    }

private:
    const Fragment* f_;
    const FragmentDistance* distance_;
};

Tree* PrintTree::make_tree(const Block* block,
        const std::string& method) const {
    Tree* tree = new Tree;
    BOOST_FOREACH (const Fragment* f, *block) {
        tree->add_node(new FragmentLeaf(f, distance_));
    }
    if (method == "upgma") {
        tree->upgma();
    } else if (method == "nj") {
        tree->neighbor_joining();
    } else {
        delete tree;
        throw Exception("Unknown tree construction method: " + method);
    }
    return tree;
}

void PrintTree::print_block(std::ostream& o, Block* block) const {
    Tree* tree = make_tree(block, opt_value("method").as<std::string>());
    o << block->name() << '\t';
    tree->print_newick(o);
    o << '\n';
    delete tree;
}

}

