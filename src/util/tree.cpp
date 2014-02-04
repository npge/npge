/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <vector>
#include <map>
#include <boost/foreach.hpp>

#include "tree.hpp"
#include "Exception.hpp"

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

float AbstractTreeNode::tree_distance_to(const AbstractTreeNode* other) const {
    typedef std::map<const AbstractTreeNode*, float> Node2Dist;
    Node2Dist node_to_dist;
    const AbstractTreeNode* node = this;
    float dist = 0.0;
    while (node) {
        node_to_dist[node] = dist;
        dist += node->length();
        node = node->parent();
    }
    node = other;
    dist = 0.0;
    while (node) {
        Node2Dist::const_iterator it = node_to_dist.find(node);
        if (it != node_to_dist.end()) {
            return it->second + dist;
        } else {
            dist += node->length();
            node = node->parent();
        }
    }
    return -1000.0;
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

float LeafNode::distance_to(const LeafNode* leaf) const {
    return distance_to_impl(leaf);
}

std::string LeafNode::name() const {
    return name_impl();
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

void Tree::add_node(AbstractTreeNode* node) {
    nodes_.insert(node);
    BranchNode* branch = dynamic_cast<BranchNode*>(node);
    if (branch) {
        add_node(branch->left());
        add_node(branch->right());
    }
}

typedef std::pair<AbstractTreeNode*, AbstractTreeNode*> Pair;
typedef std::map<Pair, float> Distances;
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

static void upgma_round(Tree* tree, Distances& distances,
        Nodes& nodes) {
    Pair min_pair;
    for (int i = 0; i < nodes.size(); i++) {
        AbstractTreeNode* node_i = nodes[i];
        for (int j = i + 1; j < nodes.size(); j++) {
            AbstractTreeNode* node_j = nodes[j];
            float distance = distances[make_pair(node_i, node_j)];
            if (min_pair == Pair() || distance < distances[min_pair]) {
                min_pair = make_pair(node_i, node_j);
            }
        }
    }
    if (!min_pair.first || !min_pair.second) {
        throw Exception("No branch for upgma round");
    }
    float min_distance = distances[min_pair];
    BranchNode* new_node = new BranchNode;
    tree->add_node(new_node);
    new_node->set_left(min_pair.first);
    new_node->set_right(min_pair.second);
    min_pair.first->set_length(min_distance / 2);
    min_pair.second->set_length(min_distance / 2);
    for (int i = 0; i < nodes.size(); i++) {
        AbstractTreeNode* node_i = nodes[i];
        if (node_i != min_pair.first && node_i != min_pair.second) {
            float d1 = distances[make_pair(node_i, min_pair.first)];
            float d2 = distances[make_pair(node_i, min_pair.second)];
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
    BOOST_FOREACH (AbstractTreeNode* node, nodes_) {
        LeafNode* leaf = dynamic_cast<LeafNode*>(node);
        if (!leaf) {
            throw Exception("Tree node is not leaf");
        }
        if (leaf->parent()) {
            throw Exception("Tree node " + leaf->name() + " has parent");
        }
        leafs.push_back(leaf);
    }
    Distances distances;
    for (int i = 0; i < leafs.size(); i++) {
        LeafNode* leaf_i = leafs[i];
        for (int j = i + 1; j < leafs.size(); j++) {
            LeafNode* leaf_j = leafs[j];
            float distance = leaf_i->distance_to(leaf_j);
            distances[make_pair(leaf_i, leaf_j)] = distance;
        }
    }
    std::vector<AbstractTreeNode*> nodes(leafs.begin(), leafs.end());
    while (!root()) {
        upgma_round(this, distances, nodes);
    }
}

}

