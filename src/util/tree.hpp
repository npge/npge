/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_TREE_HPP_
#define BR_TREE_HPP_

#include <iosfwd>
#include <set>
#include <vector>
#include <boost/utility.hpp>

namespace bloomrepeats {

class AbstractTreeNode;
class BranchNode;
class LeafNode;
class Tree;

class AbstractTreeNode : boost::noncopyable {
public:
    AbstractTreeNode();
    virtual ~AbstractTreeNode();

    AbstractTreeNode* clone() const;

    double length() const {
        return length_;
    }

    void set_length(double length) {
        length_ = length;
    }

    BranchNode* parent() const {
        return parent_;
    }

    void set_parent(BranchNode* parent) {
        parent_ = parent;
    }

    /** Distance between nodes according to tree branches lengthes.
    If nodes are unrelated, return sum of distances to roots.
    */
    double tree_distance_to(const AbstractTreeNode* other) const;

    void print_newick(std::ostream& o, bool lengthes = true) const;

    std::string newick(bool lengthes = true) const;

protected:
    virtual AbstractTreeNode* clone_impl() const = 0;
    virtual void print_newick_impl(std::ostream& o, bool lengthes) const = 0;

    double length_;
    BranchNode* parent_;
};

class BranchNode : public AbstractTreeNode {
public:
    BranchNode();

    AbstractTreeNode* left() const {
        return left_;
    }

    void set_left(AbstractTreeNode* left);

    AbstractTreeNode* right() const {
        return right_;
    }

    void set_right(AbstractTreeNode* right);

protected:
    AbstractTreeNode* clone_impl() const;
    void print_newick_impl(std::ostream& o, bool lengthes) const;

private:
    AbstractTreeNode* left_;
    AbstractTreeNode* right_;
};

class LeafNode : public AbstractTreeNode {
public:
    double distance_to(const LeafNode* leaf) const;
    std::string name() const;

protected:
    virtual double distance_to_impl(const LeafNode* leaf) const = 0;
    virtual std::string name_impl() const = 0;
    void print_newick_impl(std::ostream& o, bool lengthes) const;
};

class Tree : boost::noncopyable {
public:
    Tree();
    virtual ~Tree();

    void clear();

    AbstractTreeNode* root() const {
        return root_;
    }

    void set_root(AbstractTreeNode* node) {
        root_ = node;
    }

    const std::set<AbstractTreeNode*>& nodes() const {
        return nodes_;
    }

    std::vector<AbstractTreeNode*> orphan_nodes() const;

    Tree* clone() const;

    void print_newick(std::ostream& o, bool lengthes = true) const;

    std::string newick(bool lengthes = true) const;

    void add_node(AbstractTreeNode* node);

    void upgma();

    void neighbor_joining();

private:
    AbstractTreeNode* root_;
    std::set<AbstractTreeNode*> nodes_;
};

}

#endif

