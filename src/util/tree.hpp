/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_TREE_HPP_
#define BR_TREE_HPP_

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

    float length() const {
        return length_;
    }

    void set_length(float length) {
        length_ = length;
    }

    BranchNode* parent() const {
        return parent_;
    }

    void set_parent(BranchNode* parent) {
        parent_ = parent;
    }

    /** Distance between nodes according to tree branches lengthes.
    If nodes are unrelated, return negative number.
    */
    float tree_distance_to(const AbstractTreeNode* other) const;

protected:
    virtual AbstractTreeNode* clone_impl() const = 0;

    float length_;
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

private:
    AbstractTreeNode* left_;
    AbstractTreeNode* right_;
};

class LeafNode : public AbstractTreeNode {
public:
    float distance_to(const LeafNode* leaf) const;
    std::string name() const;

protected:
    virtual float distance_to_impl(const LeafNode* leaf) const = 0;
    virtual std::string name_impl() const = 0;
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

    void add_node(AbstractTreeNode* node);

    void upgma();

private:
    AbstractTreeNode* root_;
    std::set<AbstractTreeNode*> nodes_;
};

}

#endif

