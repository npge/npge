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

class TreeNode;
class LeafNode;

typedef std::vector<TreeNode*> Nodes;
typedef std::vector<LeafNode*> Leafs;

class TreeNode : boost::noncopyable {
public:
    TreeNode();

    virtual ~TreeNode();

    void clear();

    TreeNode* clone() const;

    double length() const {
        return length_;
    }

    void set_length(double length) {
        length_ = length;
    }

    TreeNode* parent() const {
        return parent_;
    }

    void set_parent(TreeNode* parent);

    const Nodes& children() const {
        return children_;
    }

    void all_descendants(Nodes& result) const;

    void all_leafs(Leafs& result) const;

    bool has_child(TreeNode* child) const;

    void add_child(TreeNode* child);

    void delete_child(TreeNode* child);

    void detach_child(TreeNode* child);

    void detach();

    /** Distance between nodes according to tree branches lengthes.
    If nodes are unrelated, return negative number.
    */
    double tree_distance_to(const TreeNode* other) const;

    void print_newick(std::ostream& o, bool lengthes = true) const;

    std::string newick(bool lengthes = true) const;

    void upgma();

    void neighbor_joining();

protected:
    virtual void print_newick_impl(std::ostream& o, bool lengthes) const;
    virtual TreeNode* clone_impl() const;

private:
    double length_;
    TreeNode* parent_;
    Nodes children_;
};

class LeafNode : public TreeNode {
public:
    double distance_to(const LeafNode* leaf) const;
    std::string name() const;

protected:
    virtual double distance_to_impl(const LeafNode* leaf) const = 0;
    virtual std::string name_impl() const = 0;
    void print_newick_impl(std::ostream& o, bool lengthes) const;
};

}

#endif

