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

namespace bloomrepeats {

TreeNode::TreeNode():
    length_(0.0), parent_(0)
{ }

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
    return new_node;
}

void TreeNode::set_parent(TreeNode* parent) {
    parent->add_child(this);
}

void TreeNode::all_descendants(Nodes& result) const {
    BOOST_FOREACH (TreeNode* child, children()) {
        result.push_back(child);
        child->all_descendants(result);
    }
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

void TreeNode::print_newick(std::ostream& o, bool lengthes) const {
    print_newick_impl(o, lengthes);
    o << ';';
}

std::string TreeNode::newick(bool lengthes) const {
    std::stringstream result;
    print_newick(result, lengthes);
    return result.str();
}

void TreeNode::print_newick_impl(std::ostream& o, bool lengthes) const {
    o << '(';
    bool first = true;
    BOOST_FOREACH (TreeNode* node, children()) {
        if (!first) {
            o << ',';
        } else {
            first = false;
        }
        node->print_newick_impl(o, lengthes);
    }
    o << ')';
    if (lengthes && parent()) {
        o << ':' << length();
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

void LeafNode::print_newick_impl(std::ostream& o, bool lengthes) const {
    o << name();
    if (lengthes && parent()) {
        o << ':' << length();
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
        a_(a), b_(b)
    { }

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
            IfRemove(min_pair.first, min_pair.second)), nodes.end());
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
    BOOST_ASSERT(nodes.size() == 1);
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
            IfRemove(min_pair.first, min_pair.second)), nodes.end());
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
    BOOST_ASSERT(nodes.size() == 3);
    Pair pair01 = make_pair(nodes[0], nodes[1]);
    double l0 = distance_to_first(pair01, distances, nodes);
    double l1 = distances[pair01] - l0;
    double l2 = distance_to_pair(pair01, distances, nodes[2]);
    pair01.first->set_length(l0);
    pair01.second->set_length(l1);
    nodes[2]->set_length(l2);
}

}

