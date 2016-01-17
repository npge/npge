/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_TREE_HPP_
#define NPGE_TREE_HPP_

#include <iosfwd>
#include <map>
#include <vector>
#include <boost/utility.hpp>

namespace npge {

class TreeNode;
class LeafNode;

typedef std::vector<TreeNode*> Nodes;
typedef std::vector<LeafNode*> Leafs;
typedef std::map<std::string, double> BranchTable;
// branch is a sequence of 0 and 1 (first is 0)
// 0 and 1 mark two sets of leafs

class TreeNode : boost::noncopyable {
public:
    enum ShowBootstrap {
        NO_BOOTSTRAP = 0,
        BOOTSTRAP_IN_BRACES = 1,
        BOOTSTRAP_BEFORE_LENGTH = 2
    };

    TreeNode();

    virtual ~TreeNode();

    void clear();

    TreeNode* clone() const;

    TreeNode* clone_with_pseudo_leafs() const;

    double length() const {
        return length_;
    }

    void set_length(double length) {
        length_ = length;
    }

    double bootstrap() const {
        return bootstrap_;
    }

    void set_bootstrap(double bootstrap) {
        bootstrap_ = bootstrap;
    }

    TreeNode* parent() const {
        return parent_;
    }

    void set_parent(TreeNode* parent);

    const Nodes& children() const {
        return children_;
    }

    void all_nodes(Nodes& result) const;

    void all_descendants(Nodes& result) const;

    void all_leafs_and_this(Leafs& result) const;

    void all_leafs(Leafs& result) const;

    void all_end_nodes(Nodes& result) const;

    bool has_child(TreeNode* child) const;

    void add_child(TreeNode* child);

    void delete_child(TreeNode* child);

    void detach_child(TreeNode* child);

    void detach();

    /** Distance between nodes according to tree branches lengthes.
    If nodes are unrelated, return negative number.
    */
    double tree_distance_to(const TreeNode* other) const;

    void print_newick(std::ostream& o,
                      bool lengthes = true,
                      ShowBootstrap sbs = NO_BOOTSTRAP) const;

    std::string newick(bool lengthes = true,
                       ShowBootstrap sbs = NO_BOOTSTRAP) const;

    void upgma();

    void neighbor_joining();

    void branch_table(BranchTable& table, const Leafs& leafs,
                      double weight) const;

    /** Build tree from branch information.
    Conflicting branches with lower weight  are discarded.
    */
    void from_branches(const BranchTable& table,
                       const Leafs& leafs);

    /** Return string or '0' and '1' encoding the branch.
    \param leafs All leafs.
    \param leafs Subset of leafs in branch.
    Each node divides complete leafs set (argument leafs)
    to two subsets: children and any other.
    These subsets are encoded as '0' and '1' in such a way
    that first char of resulting string is '0'.
    Order of chars in string and its length corresponds
    to argument leafs.
    */
    static std::string branch_str_encode(const Leafs& leafs,
                                         const Leafs& sub_leafs);

    /** Return string or '0' and '1' encoding the branch */
    std::string branch_str_encode(const Leafs& leafs) const;

    /** Decode branch string to two subsets of leafs */
    static void branch_str_decode(const Leafs& leafs,
                                  const std::string& branch_str,
                                  Leafs& sub_leafs_0, Leafs& sub_leafs_1);

    /** Return string with two sets.
    Example: {a, b, c} vs {e, f} */
    static std::string branch_as_sets(const Leafs& leafs,
                                      const std::string& branch_str);

    /** Returns if child branch is inside parent branch.
    If parent[i] == X for each i, then returns X.
    Returns '0' if child[i] == const X for each i
    such that parent[i] == '1'.
    Returns '1' if child[i] == const X for each i
    such that parent[i] == '0'.
    Returns '\0' otherwise, which means that branches
    are incompatible.

    Branches must not be equal!
    */
    static char branch_inside(const std::string& child,
                              const std::string& parent);

    /** Return if tree possible with both branches */
    static bool branches_compatible(const std::string& b1,
                                    const std::string& b2);

protected:
    virtual void print_newick_impl(std::ostream& o,
                                   bool lengthes,
                                   ShowBootstrap sbs) const;
    virtual TreeNode* clone_impl() const;

private:
    TreeNode* parent_;
    Nodes children_;
    double length_;
    double bootstrap_;
};

class LeafNode : public TreeNode {
public:
    double distance_to(const LeafNode* leaf) const;
    std::string name() const;

protected:
    virtual double distance_to_impl(const LeafNode* leaf) const = 0;
    virtual std::string name_impl() const = 0;
    void print_newick_impl(std::ostream& o,
                           bool lengthes,
                           ShowBootstrap sbs) const;
};

}

#endif

