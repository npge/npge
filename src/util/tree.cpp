/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <vector>
#include <sstream>
#include <map>
#include <boost/foreach.hpp>

#include "tree.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"

// TODO inherit Tree from AbstractTreeNode
// TODO multichild nodes (as Tree); use in end of neighbor_joining

namespace bloomrepeats {

AbstractTreeNode::AbstractTreeNode():
    length_(0.0), parent_(0)
{ }

AbstractTreeNode::~AbstractTreeNode()
{ }

AbstractTreeNode* AbstractTreeNode::clone() const {
    AbstractTreeNode* new_node = clone_impl();
    new_node->set_length(length());
    return new_node;
}

double AbstractTreeNode::tree_distance_to(const AbstractTreeNode* other) const {
    typedef std::map<const AbstractTreeNode*, double> Node2Dist;
    Node2Dist node_to_dist;
    const AbstractTreeNode* node = this;
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
    return dist + dist2;
}

void AbstractTreeNode::print_newick(std::ostream& o, bool lengthes) const {
    print_newick_impl(o, lengthes);
    if (lengthes && parent()) {
        o << ':' << length();
    }
}

std::string AbstractTreeNode::newick(bool lengthes) const {
    std::stringstream result;
    print_newick(result, lengthes);
    return result.str();
}

BranchNode::BranchNode():
    left_(0), right_(0)
{ }

void BranchNode::set_left(AbstractTreeNode* left) {
    left_ = left;
    if (left_) {
        left_->set_parent(this);
    }
}

void BranchNode::set_right(AbstractTreeNode* right) {
    right_ = right;
    if (right_) {
        right_->set_parent(this);
    }
}

AbstractTreeNode* BranchNode::clone_impl() const {
    BranchNode* new_node = new BranchNode;
    if (left()) {
        new_node->set_left(left()->clone());
    }
    if (right()) {
        new_node->set_right(right()->clone());
    }
    return new_node;
}

void BranchNode::print_newick_impl(std::ostream& o, bool lengthes) const {
    o << '(';
    if (left()) {
        left()->print_newick(o, lengthes);
    }
    o << ',';
    if (right()) {
        right()->print_newick(o, lengthes);
    }
    o << ')';
}

double LeafNode::distance_to(const LeafNode* leaf) const {
    return distance_to_impl(leaf);
}

std::string LeafNode::name() const {
    return name_impl();
}

void LeafNode::print_newick_impl(std::ostream& o, bool lengthes) const {
    o << name();
}

Tree::Tree():
    root_(0)
{ }

Tree::~Tree() {
    clear();
}

void Tree::clear() {
    BOOST_FOREACH (AbstractTreeNode* node, nodes_) {
        delete node;
    }
    set_root(0);
    nodes_.clear();
}

std::vector<AbstractTreeNode*> Tree::orphan_nodes() const {
    std::vector<AbstractTreeNode*> result;
    BOOST_FOREACH (AbstractTreeNode* node, nodes()) {
        if (!node->parent()) {
            result.push_back(node);
        }
    }
    return result;
}

Tree* Tree::clone() const {
    Tree* new_tree = new Tree;
    BOOST_FOREACH (AbstractTreeNode* node, orphan_nodes()) {
        new_tree->add_node(node->clone());
    }
    return new_tree;
}

void Tree::print_newick(std::ostream& o, bool lengthes) const {
    std::vector<AbstractTreeNode*> orphans = orphan_nodes();
    if (orphans.size() >= 2) {
        o << '(';
    }
    bool first = true;
    BOOST_FOREACH (AbstractTreeNode* node, orphans) {
        if (!first) {
            o << ',';
        } else {
            first = false;
        }
        node->print_newick(o, lengthes);
    }
    if (orphans.size() >= 2) {
        o << ')';
    }
    o << ';';
}

std::string Tree::newick(bool lengthes) const {
    std::stringstream result;
    print_newick(result, lengthes);
    return result.str();
}

void Tree::add_node(AbstractTreeNode* node) {
    BOOST_ASSERT(node);
    nodes_.insert(node);
    BranchNode* branch = dynamic_cast<BranchNode*>(node);
    if (branch) {
        if (branch->left()) {
            add_node(branch->left());
        }
        if (branch->right()) {
            add_node(branch->right());
        }
    }
}

typedef std::pair<AbstractTreeNode*, AbstractTreeNode*> Pair;
typedef std::map<Pair, double> Distances;
typedef std::vector<AbstractTreeNode*> Nodes;

Pair make_pair(AbstractTreeNode* a, AbstractTreeNode* b) {
    if (a < b) {
        return Pair(a, b);
    } else {
        return Pair(b, a);
    }
}

struct IfRemove {
    IfRemove(AbstractTreeNode* a, AbstractTreeNode* b):
        a_(a), b_(b)
    { }

    bool operator()(AbstractTreeNode* node) {
        return node == a_ || node == b_;
    }

private:
    AbstractTreeNode* a_;
    AbstractTreeNode* b_;
};

static void find_leafs_and_distances(Tree* tree,
        std::vector<LeafNode*>& leafs,
        Distances& distances) {
    BOOST_FOREACH (AbstractTreeNode* node, tree->nodes()) {
        LeafNode* leaf = dynamic_cast<LeafNode*>(node);
        if (!leaf) {
            throw Exception("Tree node is not leaf");
        }
        if (leaf->parent()) {
            throw Exception("Tree node " + leaf->name() + " has parent");
        }
        leafs.push_back(leaf);
    }
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
        AbstractTreeNode* node_i = nodes[i];
        for (int j = i + 1; j < nodes.size(); j++) {
            AbstractTreeNode* node_j = nodes[j];
            double distance = distances[make_pair(node_i, node_j)];
            if (min_pair == Pair() || distance < distances[min_pair]) {
                min_pair = make_pair(node_i, node_j);
            }
        }
    }
    return min_pair;
}

static void upgma_round(Tree* tree, Distances& distances,
        Nodes& nodes) {
    Pair min_pair = find_min_pair(distances, nodes);
    if (!min_pair.first || !min_pair.second) {
        throw Exception("No branch for upgma round");
    }
    double min_distance = distances[min_pair];
    BranchNode* new_node = new BranchNode;
    tree->add_node(new_node);
    new_node->set_left(min_pair.first);
    new_node->set_right(min_pair.second);
    min_pair.first->set_length(min_distance / 2.0);
    min_pair.second->set_length(min_distance / 2.0);
    for (int i = 0; i < nodes.size(); i++) {
        AbstractTreeNode* node_i = nodes[i];
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
            IfRemove(min_pair.first, min_pair.second)), nodes.end());
    nodes.push_back(new_node);
    if (distances.empty()) {
        tree->set_root(new_node);
    }
}

void Tree::upgma() {
    std::vector<LeafNode*> leafs;
    Distances distances;
    find_leafs_and_distances(this, leafs, distances);
    std::vector<AbstractTreeNode*> nodes(leafs.begin(), leafs.end());
    while (!root()) {
        upgma_round(this, distances, nodes);
    }
}

static void calculate_q(Distances& Q, Distances& distances,
        Nodes& nodes) {
    for (int i = 0; i < nodes.size(); i++) {
        AbstractTreeNode* node_i = nodes[i];
        for (int j = i + 1; j < nodes.size(); j++) {
            AbstractTreeNode* node_j = nodes[j];
            double distance = (nodes.size() - 2.0) *
                distances[make_pair(node_i, node_j)];
            for (int k = 0; k < nodes.size(); k++) {
                AbstractTreeNode* node_k = nodes[k];
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
    AbstractTreeNode* node_i = min_pair.first;
    AbstractTreeNode* node_j = min_pair.second;
    for (int k = 0; k < nodes.size(); k++) {
        AbstractTreeNode* node_k = nodes[k];
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
        Distances& distances, AbstractTreeNode* node_k) {
    double min_distance = distances[min_pair];
    double dist_f = distances[make_pair(min_pair.first, node_k)];
    double dist_s = distances[make_pair(min_pair.second, node_k)];
    return 0.5 * (dist_f + dist_s - min_distance);
}

static void neighbor_joining_round(Tree* tree, Distances& distances,
        Nodes& nodes) {
    Distances Q;
    calculate_q(Q, distances, nodes);
    Pair min_pair = find_min_pair(Q, nodes);
    if (!min_pair.first || !min_pair.second) {
        throw Exception("No min element of Q for neighbor joining");
    }
    BranchNode* new_node = new BranchNode;
    tree->add_node(new_node);
    new_node->set_left(min_pair.first);
    new_node->set_right(min_pair.second);
    double min_distance = distances[min_pair];
    double distance_to_left = distance_to_first(min_pair, distances, nodes);
    double distance_to_right = min_distance - distance_to_left;
    min_pair.first->set_length(distance_to_left);
    min_pair.second->set_length(distance_to_right);
    for (int k = 0; k < nodes.size(); k++) {
        AbstractTreeNode* node_k = nodes[k];
        if (node_k != min_pair.first && node_k != min_pair.second) {
            double d = distance_to_pair(min_pair, distances, node_k);
            distances[make_pair(new_node, node_k)] = d;
            distances.erase(make_pair(min_pair.first, node_k));
            distances.erase(make_pair(min_pair.second, node_k));
        }
    }
    distances.erase(min_pair);
    nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
            IfRemove(min_pair.first, min_pair.second)), nodes.end());
    nodes.push_back(new_node);
}

void Tree::neighbor_joining() {
    std::vector<LeafNode*> leafs;
    Distances distances;
    find_leafs_and_distances(this, leafs, distances);
    std::vector<AbstractTreeNode*> nodes(leafs.begin(), leafs.end());
    for (int round = 0; round < leafs.size() - 3; round++) {
        neighbor_joining_round(this, distances, nodes);
    }
    BOOST_ASSERT(nodes.size() == 3);
    Pair pair01 = make_pair(nodes[0], nodes[1]);
    double l0 = distance_to_first(pair01, distances, nodes);
    double l1 = distances[pair01] - l0;
    double l2 = distance_to_pair(pair01, distances, nodes[2]);
    nodes[0]->set_length(l0);
    nodes[1]->set_length(l1);
    nodes[2]->set_length(l2);
}

}

