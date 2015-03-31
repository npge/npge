/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <vector>
#include <sstream>
#include <set>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>

#include "tree.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"

namespace npge {

TreeNode::TreeNode():
    parent_(0), length_(0.0), bootstrap_(0.0) {
}

TreeNode::~TreeNode() {
    clear();
    detach();
}

void TreeNode::clear() {
    BOOST_FOREACH (TreeNode* node, children()) {
        node->parent_ = 0;
        delete node;
    }
    children_.clear();
}

TreeNode* TreeNode::clone() const {
    TreeNode* new_node = clone_impl();
    new_node->set_length(length());
    new_node->set_bootstrap(bootstrap());
    return new_node;
}

class DotNode : public LeafNode {
public:
    DotNode(const LeafNode* leaf) {
        name_ = "." + leaf->name();
    }

protected:
    double distance_to_impl(const LeafNode* leaf) const {
        // useless
        return 0.0;
    }

    std::string name_impl() const {
        return name_;
    }

private:
    std::string name_;
};

TreeNode* TreeNode::clone_with_pseudo_leafs() const {
    TreeNode* new_node = new TreeNode;
    const LeafNode* leaf = dynamic_cast<const LeafNode*>(this);
    if (leaf) {
        TreeNode* c = leaf->clone();
        c->set_length(0.0);
        c->set_bootstrap(0.0);
        new_node->add_child(c);
        new_node->add_child(new DotNode(leaf));
    } else {
        BOOST_FOREACH (TreeNode* child, children()) {
            new_node->add_child(child->clone_with_pseudo_leafs());
        }
    }
    new_node->set_length(length());
    new_node->set_bootstrap(bootstrap());
    return new_node;
}

void TreeNode::set_parent(TreeNode* parent) {
    parent->add_child(this);
}

void TreeNode::all_nodes(Nodes& result) const {
    result.push_back(const_cast<TreeNode*>(this));
    all_descendants(result);
}

void TreeNode::all_descendants(Nodes& result) const {
    BOOST_FOREACH (TreeNode* child, children()) {
        result.push_back(child);
        child->all_descendants(result);
    }
}

void TreeNode::all_leafs_and_this(Leafs& result) const {
    const LeafNode* leaf = dynamic_cast<const LeafNode*>(this);
    if (leaf) {
        result.push_back(const_cast<LeafNode*>(leaf));
    }
    all_leafs(result);
}

void TreeNode::all_leafs(Leafs& result) const {
    BOOST_FOREACH (TreeNode* child, children()) {
        LeafNode* leaf = dynamic_cast<LeafNode*>(child);
        if (leaf) {
            result.push_back(leaf);
        }
        child->all_leafs(result);
    }
}

void TreeNode::all_end_nodes(Nodes& result) const {
    if (children().empty()) {
        result.push_back(const_cast<TreeNode*>(this));
    } else {
        BOOST_FOREACH (TreeNode* child, children()) {
            child->all_end_nodes(result);
        }
    }
}

bool TreeNode::has_child(TreeNode* child) const {
    BOOST_FOREACH (TreeNode* node, children()) {
        if (node == child) {
            return true;
        }
    }
    return false;
}

void TreeNode::add_child(TreeNode* child) {
    if (child->parent()) {
        child->detach();
    }
    if (!has_child(child)) {
        children_.push_back(child);
    }
    child->parent_ = this;
}

void TreeNode::delete_child(TreeNode* child) {
    detach_child(child);
    delete child;
}

void TreeNode::detach_child(TreeNode* child) {
    children_.erase(std::remove(children_.begin(), children_.end(),
                                child), children_.end());
    child->parent_ = 0;
}

void TreeNode::detach() {
    if (parent()) {
        parent()->detach_child(this);
    }
}

double TreeNode::tree_distance_to(const TreeNode* other) const {
    typedef std::map<const TreeNode*, double> Node2Dist;
    Node2Dist node_to_dist;
    const TreeNode* node = this;
    double dist = 0.0;
    while (node) {
        node_to_dist[node] = dist;
        dist += node->length();
        node = node->parent();
    }
    node = other;
    double dist2 = 0.0;
    while (node) {
        Node2Dist::const_iterator it = node_to_dist.find(node);
        if (it != node_to_dist.end()) {
            return it->second + dist2;
        } else {
            dist2 += node->length();
            node = node->parent();
        }
    }
    return -1000.0;
}

void TreeNode::print_newick(std::ostream& o, bool lengthes,
                            ShowBootstrap sbs) const {
    print_newick_impl(o, lengthes, sbs);
    o << ';';
}

std::string TreeNode::newick(bool lengthes,
                             ShowBootstrap sbs) const {
    std::stringstream result;
    print_newick(result, lengthes, sbs);
    return result.str();
}

void TreeNode::print_newick_impl(std::ostream& o, bool lengthes,
                                 ShowBootstrap sbs) const {
    int ancestors_number = 0;
    TreeNode* p = parent();
    while (p) {
        ancestors_number += 1;
        p = p->parent();
    }
    std::string indent(ancestors_number * 2, ' ');
    o << '(';
    bool first = true;
    BOOST_FOREACH (TreeNode* node, children()) {
        if (!first) {
            o << ',';
        } else {
            first = false;
        }
        o << '\n' << indent << "  ";
        node->print_newick_impl(o, lengthes, sbs);
    }
    o << "\n" << indent << ')';
    if (lengthes && parent()) {
        if (sbs == BOOTSTRAP_BEFORE_LENGTH) {
            o << bootstrap();
        }
        o << ':' << length();
    }
    if (sbs == BOOTSTRAP_IN_BRACES) {
        o << '[' << bootstrap() << ']';
    }
}

TreeNode* TreeNode::clone_impl() const {
    TreeNode* new_node = new TreeNode;
    BOOST_FOREACH (TreeNode* child, children()) {
        new_node->add_child(child->clone());
    }
    return new_node;
}

double LeafNode::distance_to(const LeafNode* leaf) const {
    return distance_to_impl(leaf);
}

std::string LeafNode::name() const {
    return name_impl();
}

void LeafNode::print_newick_impl(std::ostream& o, bool lengthes,
                                 ShowBootstrap sbs) const {
    o << name();
    if (lengthes && parent()) {
        o << ':' << length();
    }
    if (sbs == BOOTSTRAP_IN_BRACES) {
        o << '[' << bootstrap() << ']';
    }
}

typedef std::pair<TreeNode*, TreeNode*> Pair;
typedef std::map<Pair, double> Distances;

Pair make_pair(TreeNode* a, TreeNode* b) {
    if (a < b) {
        return Pair(a, b);
    } else {
        return Pair(b, a);
    }
}

struct IfRemove {
    IfRemove(TreeNode* a, TreeNode* b):
        a_(a), b_(b) {
    }

    bool operator()(TreeNode* node) {
        return node == a_ || node == b_;
    }

private:
    TreeNode* a_;
    TreeNode* b_;
};

static void build_distances(Distances& distances, Leafs& leafs) {
    for (int i = 0; i < leafs.size(); i++) {
        LeafNode* leaf_i = leafs[i];
        for (int j = i + 1; j < leafs.size(); j++) {
            LeafNode* leaf_j = leafs[j];
            double distance = leaf_i->distance_to(leaf_j);
            distances[make_pair(leaf_i, leaf_j)] = distance;
        }
    }
}

static Pair find_min_pair(Distances& distances, const Nodes& nodes) {
    Pair min_pair;
    for (int i = 0; i < nodes.size(); i++) {
        TreeNode* node_i = nodes[i];
        for (int j = i + 1; j < nodes.size(); j++) {
            TreeNode* node_j = nodes[j];
            double distance = distances[make_pair(node_i, node_j)];
            if (min_pair == Pair() || distance < distances[min_pair]) {
                min_pair = make_pair(node_i, node_j);
            }
        }
    }
    return min_pair;
}

static void upgma_round(TreeNode* tree, Distances& distances,
                        Nodes& nodes) {
    Pair min_pair = find_min_pair(distances, nodes);
    if (!min_pair.first || !min_pair.second) {
        throw Exception("No branch for upgma round");
    }
    double min_distance = distances[min_pair];
    TreeNode* new_node = new TreeNode;
    tree->add_child(new_node);
    new_node->add_child(min_pair.first);
    new_node->add_child(min_pair.second);
    min_pair.first->set_length(min_distance / 2.0);
    min_pair.second->set_length(min_distance / 2.0);
    for (int i = 0; i < nodes.size(); i++) {
        TreeNode* node_i = nodes[i];
        if (node_i != min_pair.first && node_i != min_pair.second) {
            double d1 = distances[make_pair(node_i, min_pair.first)];
            double d2 = distances[make_pair(node_i, min_pair.second)];
            distances[make_pair(node_i, new_node)] = 0.5 * d1 + 0.5 * d2;
            distances.erase(make_pair(node_i, min_pair.first));
            distances.erase(make_pair(node_i, min_pair.second));
        }
    }
    distances.erase(min_pair);
    nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
                               IfRemove(min_pair.first, min_pair.second)),
                nodes.end());
    nodes.push_back(new_node);
}

void TreeNode::upgma() {
    Leafs leafs;
    all_leafs(leafs);
    Distances distances;
    build_distances(distances, leafs);
    BOOST_FOREACH (LeafNode* leaf, leafs) {
        leaf->detach();
    }
    clear();
    BOOST_FOREACH (LeafNode* leaf, leafs) {
        add_child(leaf);
    }
    Nodes nodes(leafs.begin(), leafs.end());
    for (int round = 0; round < int(leafs.size()) - 1; round++) {
        upgma_round(this, distances, nodes);
    }
    ASSERT_EQ(nodes.size(), 1);
}

static void calculate_q(Distances& Q, Distances& distances, Nodes& nodes) {
    for (int i = 0; i < nodes.size(); i++) {
        TreeNode* node_i = nodes[i];
        for (int j = i + 1; j < nodes.size(); j++) {
            TreeNode* node_j = nodes[j];
            double distance = (nodes.size() - 2.0) *
                              distances[make_pair(node_i, node_j)];
            for (int k = 0; k < nodes.size(); k++) {
                TreeNode* node_k = nodes[k];
                distance -= distances[make_pair(node_i, node_k)];
                distance -= distances[make_pair(node_j, node_k)];
            }
            Q[make_pair(node_i, node_j)] = distance;
        }
    }
}

static double distance_to_first(const Pair& min_pair,
                                Distances& distances, Nodes& nodes) {
    double min_distance = distances[min_pair];
    double s = 0;
    TreeNode* node_i = min_pair.first;
    TreeNode* node_j = min_pair.second;
    for (int k = 0; k < nodes.size(); k++) {
        TreeNode* node_k = nodes[k];
        if (node_k != node_i && node_k != node_j) {
            s += distances[make_pair(node_i, node_k)];
            s -= distances[make_pair(node_j, node_k)];
        }
    }
    double dist;
    if (nodes.size() > 2) {
        dist = 0.5 * min_distance + 0.5 * s / (nodes.size() - 2);
    } else {
        dist = 0.5 * min_distance;
    }
    if (dist < 0.0) {
        dist = 0.0;
    }
    if (dist > min_distance) {
        dist = min_distance;
    }
    return dist;
}

static double distance_to_pair(const Pair& min_pair,
                               Distances& distances, TreeNode* node_k) {
    double min_distance = distances[min_pair];
    double dist_f = distances[make_pair(min_pair.first, node_k)];
    double dist_s = distances[make_pair(min_pair.second, node_k)];
    return 0.5 * (dist_f + dist_s - min_distance);
}

static void neighbor_joining_round(TreeNode* tree, Distances& distances,
                                   Nodes& nodes) {
    Distances Q;
    calculate_q(Q, distances, nodes);
    Pair min_pair = find_min_pair(Q, nodes);
    if (!min_pair.first || !min_pair.second) {
        throw Exception("No min element of Q for neighbor joining");
    }
    TreeNode* new_node = new TreeNode;
    tree->add_child(new_node);
    new_node->add_child(min_pair.first);
    new_node->add_child(min_pair.second);
    double min_distance = distances[min_pair];
    double distance_to_left = distance_to_first(min_pair, distances, nodes);
    double distance_to_right = min_distance - distance_to_left;
    min_pair.first->set_length(distance_to_left);
    min_pair.second->set_length(distance_to_right);
    for (int k = 0; k < nodes.size(); k++) {
        TreeNode* node_k = nodes[k];
        if (node_k != min_pair.first && node_k != min_pair.second) {
            double d = distance_to_pair(min_pair, distances, node_k);
            distances[make_pair(new_node, node_k)] = d;
            distances.erase(make_pair(min_pair.first, node_k));
            distances.erase(make_pair(min_pair.second, node_k));
        }
    }
    distances.erase(min_pair);
    nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
                               IfRemove(min_pair.first, min_pair.second)),
                nodes.end());
    nodes.push_back(new_node);
}

void TreeNode::neighbor_joining() {
    Leafs leafs;
    all_leafs(leafs);
    Distances distances;
    build_distances(distances, leafs);
    BOOST_FOREACH (LeafNode* leaf, leafs) {
        leaf->detach();
    }
    clear();
    BOOST_FOREACH (LeafNode* leaf, leafs) {
        add_child(leaf);
    }
    Nodes nodes(leafs.begin(), leafs.end());
    if (nodes.size() == 0) {
        return;
    } else if (nodes.size() == 1) {
        return;
    } else if (nodes.size() == 2) {
        double d = distances[make_pair(nodes[0], nodes[1])];
        nodes[0]->set_length(d / 2.0);
        nodes[1]->set_length(d / 2.0);
        return;
    }
    for (int round = 0; round < int(leafs.size()) - 3; round++) {
        neighbor_joining_round(this, distances, nodes);
    }
    ASSERT_EQ(nodes.size(), 3);
    Pair pair01 = make_pair(nodes[0], nodes[1]);
    double l0 = distance_to_first(pair01, distances, nodes);
    double l1 = distances[pair01] - l0;
    double l2 = distance_to_pair(pair01, distances, nodes[2]);
    pair01.first->set_length(l0);
    pair01.second->set_length(l1);
    nodes[2]->set_length(l2);
}

void TreeNode::branch_table(BranchTable& table, const Leafs& leafs,
                            double weight) const {
    Nodes nodes;
    all_descendants(nodes);
    BOOST_FOREACH (TreeNode* branch, nodes) {
        Leafs sub_leafs;
        branch->all_leafs(sub_leafs);
        if (sub_leafs.size() >= 2 && leafs.size() - sub_leafs.size() >= 2) {
            std::string branch_str = branch_str_encode(leafs, sub_leafs);
            double branch_weight = branch->length() * weight;
            table[branch_str] += branch_weight;
        }
    }
}

static std::string make_subbranch(const std::string& child,
        const std::string& parent, char parent_char) {
    ASSERT_EQ(child.size(), parent.size());
    ASSERT_NE(child, parent);
    ASSERT_TRUE(parent_char == '0' || parent_char == '1');
    std::string result;
    for (int i = 0; i < child.size(); i++) {
        if (parent[i] == parent_char) {
            result.push_back(child[i]);
        }
    }
    ASSERT_GT(result.size(), 0);
    if (result[0] == '1') {
        BOOST_FOREACH (char& c, result) {
            c = (c == '0') ? '1' : '0';
        }
    }
    return result;
}

void TreeNode::from_branches(const BranchTable& table,
                             const Leafs& leafs) {
    if (table.empty()) {
        BOOST_FOREACH (LeafNode* leaf, leafs) {
            add_child(leaf);
        }
        return;
    }
    std::string best_branch;
    double best_weight = -1;
    BOOST_FOREACH (const BranchTable::value_type& v, table) {
        const std::string& branch = v.first;
        const double& weight = v.second;
        ASSERT_GT(branch.size(), 0);
        ASSERT_EQ(branch[0], '0');
        if (weight > best_weight) {
            best_branch = branch;
            best_weight = weight;
        }
    }
    Leafs sub_leafs_0, sub_leafs_1;
    branch_str_decode(leafs, best_branch,
                      sub_leafs_0, sub_leafs_1);
    ASSERT_EQ(sub_leafs_0.size() + sub_leafs_1.size(),
              leafs.size());
    ASSERT_GT(sub_leafs_0.size(), 0);
    ASSERT_GT(sub_leafs_1.size(), 0);
    // make new branch tables
    BranchTable compatible_branches0, compatible_branches1;
    BOOST_FOREACH (const BranchTable::value_type& v, table) {
        const std::string& branch = v.first;
        const double& weight = v.second;
        if (branch != best_branch) {
            char comp = branch_inside(branch, best_branch);
            if (comp) {
                BranchTable& cb = (comp == '0') ?
                    (compatible_branches0) :
                    (compatible_branches1);
                std::string child = make_subbranch(branch,
                                    best_branch, comp);
                cb[child] = weight;
            }
        }
    }
    if (sub_leafs_0.size() == 1) {
        add_child(sub_leafs_0[0]);
    } else {
        TreeNode* node0 = new TreeNode;
        add_child(node0);
        node0->from_branches(compatible_branches0, sub_leafs_0);
    }
    if (sub_leafs_1.size() == 1) {
        add_child(sub_leafs_1[0]);
    } else {
        TreeNode* node1 = new TreeNode;
        add_child(node1);
        node1->from_branches(compatible_branches1, sub_leafs_1);
    }
}

std::string TreeNode::branch_str_encode(const Leafs& leafs,
                                        const Leafs& sub_leafs) {
    std::set<LeafNode*> sub_leafs_set(sub_leafs.begin(),
                                      sub_leafs.end());
    // ol is 0 or 1
    unsigned ol = (sub_leafs_set.find(leafs[0]) != sub_leafs_set.end());
    std::string result;
    result.reserve(leafs.size());
    BOOST_FOREACH (LeafNode* leaf, leafs) {
        // l is 0 or 1
        unsigned l = (sub_leafs_set.find(leaf) != sub_leafs_set.end());
        unsigned lo = 1 & (ol ^ l);
        result += (lo ? '1' : '0');
    }
    return result;
}

std::string TreeNode::branch_str_encode(const Leafs& leafs) const {
    Leafs sub_leafs;
    all_leafs(sub_leafs);
    return branch_str_encode(leafs, sub_leafs);
}

void TreeNode::branch_str_decode(const Leafs& leafs,
                                 const std::string& branch_str,
                                 Leafs& sub_leafs_0, Leafs& sub_leafs_1) {
    ASSERT_EQ(leafs.size(), branch_str.size());
    for (int i = 0; i < leafs.size(); i++) {
        LeafNode* leaf = leafs[i];
        if (branch_str[i] == '0') {
            sub_leafs_0.push_back(leaf);
        } else {
            sub_leafs_1.push_back(leaf);
        }
    }
}

std::string TreeNode::branch_as_sets(const Leafs& leafs,
                                     const std::string& branch_str) {
    Leafs s0, s1;
    branch_str_decode(leafs, branch_str, s0, s1);
    Strings n0, n1;
    BOOST_FOREACH (LeafNode* leaf, s0) {
        n0.push_back(leaf->name());
    }
    BOOST_FOREACH (LeafNode* leaf, s1) {
        n1.push_back(leaf->name());
    }
    using namespace boost::algorithm;
    return "{" + join(n0, ",") + "} vs {" + join(n1, ",") + "}";
}

char TreeNode::branch_inside(const std::string& child,
                             const std::string& parent) {
    ASSERT_EQ(child.size(), parent.size());
    ASSERT_NE(child, parent);
    bool seen[2][2] = {{false, false}, {false, false}};
    for (int i = 0; i < child.size(); i++) {
        int child_i = (child[i] == '0') ? 0 : 1;
        int parent_i = (parent[i] == '0') ? 0 : 1;
        seen[child_i][parent_i] = true;
    }
    if (!seen[0][0] && !seen[1][0]) {
        ASSERT_TRUE(seen[0][1] || seen[1][1]);
        return '1';
    }
    if (!seen[0][1] && !seen[1][1]) {
        ASSERT_TRUE(seen[0][0] || seen[1][0]);
        return '0';
    }
    if (seen[0][0] != seen[1][0]) {
        ASSERT_TRUE(seen[0][1] || seen[1][1]);
        return '1';
    } if (seen[0][1] != seen[1][1]) {
        ASSERT_TRUE(seen[0][0] || seen[1][0]);
        return '0';
    } else {
        return 0;
    }
}

bool TreeNode::branches_compatible(const std::string& b1,
                                   const std::string& b2) {
    return (b1 == b2) || branch_inside(b1, b2) != 0;
}

}

