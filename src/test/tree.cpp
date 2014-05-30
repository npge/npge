/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <utility>
#include <boost/test/unit_test.hpp>

#include "tree.hpp"

namespace npge {

class TestLeaf;

typedef std::pair<const LeafNode*, const LeafNode*> Pair;
typedef std::map<Pair, double> Map;

Map map;

Pair make_pair(const LeafNode* a, const LeafNode* b) {
    if (a < b) {
        return Pair(a, b);
    } else {
        return Pair(b, a);
    }
}

class TestLeaf : public LeafNode {
public:
    TestLeaf(std::string n):
        name_(n) {
    }

    double distance_to_impl(const LeafNode* leaf) const {
        return map[make_pair(this, leaf)];
    }

    std::string name_impl() const {
        return name_;
    }

    TreeNode* clone_impl() const {
        return new TestLeaf(name_);
    }

private:
    std::string name_;
};

}

bool almost_equal(double a, double b) {
    return -0.0001 < a - b && a - b < 0.0001;
}

BOOST_AUTO_TEST_CASE (tree_upgma) {
    using namespace npge;
    TreeNode tree;
    TestLeaf* a1 = new TestLeaf("a1");
    TestLeaf* a2 = new TestLeaf("a2");
    TestLeaf* a3 = new TestLeaf("a3");
    tree.add_child(a1);
    tree.add_child(a2);
    tree.add_child(a3);
    map[make_pair(a1, a2)] = 2.0;
    map[make_pair(a1, a3)] = 5.0;
    map[make_pair(a2, a3)] = 7.0;
    tree.upgma();
    BOOST_REQUIRE(tree.children().size() == 1);
    TreeNode* root = tree.children().front();
    BOOST_CHECK(almost_equal(root->length(), 0.0));
    BOOST_REQUIRE(root->children().size() == 2);
    TreeNode* left = root->children()[0];
    TreeNode* right = root->children()[1];
    BOOST_REQUIRE(left);
    BOOST_REQUIRE(right);
    BOOST_CHECK(almost_equal(left->length(), 3.0));
    BOOST_CHECK(almost_equal(right->length(), 3.0));
    BOOST_REQUIRE(dynamic_cast<TestLeaf*>(left) == a3 ||
                  dynamic_cast<TestLeaf*>(right) == a3);
    TreeNode* branch12 = (dynamic_cast<TestLeaf*>(left) == a3) ?
                         dynamic_cast<TreeNode*>(right) :
                         dynamic_cast<TreeNode*>(left);
    BOOST_REQUIRE(branch12->children().size() == 2);
    TestLeaf* br_left = dynamic_cast<TestLeaf*>(branch12->children()[0]);
    TestLeaf* br_right = dynamic_cast<TestLeaf*>(branch12->children()[1]);
    BOOST_REQUIRE(br_left);
    BOOST_REQUIRE(br_right);
    BOOST_CHECK(make_pair(br_left, br_right) == make_pair(a1, a2));
    BOOST_CHECK(almost_equal(a1->length(), 1.0));
    BOOST_CHECK(almost_equal(a2->length(), 1.0));
}

BOOST_AUTO_TEST_CASE (tree_distance) {
    using namespace npge;
    TreeNode tree;
    TestLeaf* a1 = new TestLeaf("a1");
    TestLeaf* a2 = new TestLeaf("a2");
    TestLeaf* a3 = new TestLeaf("a3");
    a1->set_length(1);
    a2->set_length(2);
    TreeNode* a12 = new TreeNode;
    a12->add_child(a1);
    a12->add_child(a2);
    a12->set_length(10);
    TreeNode* a123 = new TreeNode;
    a3->set_length(20);
    a123->add_child(a12);
    a123->add_child(a3);
    TestLeaf* a4 = new TestLeaf("a4");
    tree.add_child(a123);
    tree.add_child(a4);
    TestLeaf* a5 = new TestLeaf("a5");
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a2), 3));
    BOOST_CHECK(almost_equal(a2->tree_distance_to(a1), 3));
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a1), 0));
    BOOST_CHECK(almost_equal(a12->tree_distance_to(a12), 0));
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a12), 1));
    BOOST_CHECK(almost_equal(a2->tree_distance_to(a12), 2));
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a3), 31));
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a123), 11));
    BOOST_CHECK(almost_equal(a2->tree_distance_to(a3), 32));
    BOOST_CHECK(almost_equal(a3->tree_distance_to(a123), 20));
    BOOST_CHECK(almost_equal(a3->tree_distance_to(a12), 30));
    BOOST_CHECK(almost_equal(a3->tree_distance_to(a3), 0));
    BOOST_CHECK(almost_equal(a4->tree_distance_to(a4), 0));
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a4), 11));
    BOOST_CHECK(almost_equal(a2->tree_distance_to(a4), 12));
    BOOST_CHECK(almost_equal(a3->tree_distance_to(a4), 20));
    BOOST_CHECK(almost_equal(a123->tree_distance_to(a4), 0));
    BOOST_CHECK(a123->tree_distance_to(a5) < 0);
    BOOST_CHECK(almost_equal(a5->tree_distance_to(a5), 0));
    delete a5;
}

BOOST_AUTO_TEST_CASE (tree_nj) {
    using namespace npge;
    TreeNode tree;
    TestLeaf* a1 = new TestLeaf("a1");
    TestLeaf* a2 = new TestLeaf("a2");
    TestLeaf* a3 = new TestLeaf("a3");
    tree.add_child(a1);
    tree.add_child(a2);
    tree.add_child(a3);
    map[make_pair(a1, a2)] = 2.0;
    map[make_pair(a1, a3)] = 5.0;
    map[make_pair(a2, a3)] = 7.0;
    tree.neighbor_joining();
    BOOST_REQUIRE(tree.children().size() == 3);
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a2), 2.0));
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a3), 5.0));
    BOOST_CHECK(almost_equal(a2->tree_distance_to(a3), 7.0));
}

BOOST_AUTO_TEST_CASE (tree_nj_2) {
    using namespace npge;
    TreeNode tree;
    TestLeaf* a1 = new TestLeaf("a1");
    TestLeaf* a2 = new TestLeaf("a2");
    tree.add_child(a1);
    tree.add_child(a2);
    map[make_pair(a1, a2)] = 2.0;
    tree.neighbor_joining();
    BOOST_REQUIRE(tree.children().size() == 2);
    BOOST_CHECK(almost_equal(a1->tree_distance_to(a2), 2.0));
}

BOOST_AUTO_TEST_CASE (tree_nj_1) {
    using namespace npge;
    TreeNode tree;
    TestLeaf* a1 = new TestLeaf("a1");
    tree.add_child(a1);
    tree.neighbor_joining();
    BOOST_REQUIRE(tree.children().size() == 1);
    BOOST_CHECK(almost_equal(a1->length(), 0.0));
}

BOOST_AUTO_TEST_CASE (tree_nj_w) {
    using namespace npge;
    TreeNode tree;
    TestLeaf* a = new TestLeaf("a");
    TestLeaf* b = new TestLeaf("b");
    TestLeaf* c = new TestLeaf("c");
    TestLeaf* d = new TestLeaf("d");
    tree.add_child(a);
    tree.add_child(b);
    tree.add_child(c);
    tree.add_child(d);
    map[make_pair(a, b)] = 7.0;
    map[make_pair(a, c)] = 11.0;
    map[make_pair(a, d)] = 14.0;
    map[make_pair(b, c)] = 6.0;
    map[make_pair(b, d)] = 9.0;
    map[make_pair(c, d)] = 7.0;
    tree.neighbor_joining();
    BOOST_CHECK(almost_equal(a->tree_distance_to(b), 7.0));
    BOOST_CHECK(almost_equal(a->tree_distance_to(c), 11.0));
    BOOST_CHECK(almost_equal(a->tree_distance_to(d), 14.0));
    BOOST_CHECK(almost_equal(b->tree_distance_to(c), 6.0));
    BOOST_CHECK(almost_equal(b->tree_distance_to(d), 9.0));
    BOOST_CHECK(almost_equal(c->tree_distance_to(d), 7.0));
    BOOST_CHECK(a->parent() == b->parent());
}

BOOST_AUTO_TEST_CASE (tree_branch_str) {
    using namespace npge;
    TreeNode tree;
    TestLeaf* a1 = new TestLeaf("a1");
    TestLeaf* a2 = new TestLeaf("a2");
    a1->set_length(1);
    a2->set_length(2);
    TreeNode* a12 = new TreeNode;
    a12->add_child(a1);
    a12->add_child(a2);
    a12->set_length(10);
    TestLeaf* a3 = new TestLeaf("a3");
    TestLeaf* a4 = new TestLeaf("a4");
    a3->set_length(1);
    a4->set_length(2);
    TreeNode* a34 = new TreeNode;
    a34->add_child(a3);
    a34->add_child(a4);
    a34->set_length(20);
    tree.add_child(a12);
    tree.add_child(a34);
    Leafs leafs;
    leafs.push_back(a1);
    leafs.push_back(a2);
    leafs.push_back(a3);
    leafs.push_back(a4);
    BOOST_CHECK(a12->branch_str_encode(leafs) == "0011");
    BranchTable table;
    tree.branch_table(table, leafs, 1.0);
    BOOST_REQUIRE(table.size() == 1);
    BOOST_CHECK(table.begin()->first == "0011");
    BOOST_CHECK(almost_equal(table.begin()->second, 10.0 + 20.0));
    Leafs l0, l1;
    TreeNode::branch_str_decode(leafs, "0011", l0, l1);
    BOOST_REQUIRE(l0.size() == 2);
    BOOST_REQUIRE(l1.size() == 2);
    BOOST_CHECK(l0[0] == a1);
    BOOST_CHECK(l0[1] == a2);
    BOOST_CHECK(l1[0] == a3);
    BOOST_CHECK(l1[1] == a4);
}

BOOST_AUTO_TEST_CASE (tree_branch_str_compatible) {
    using namespace npge;
    BOOST_CHECK(TreeNode::branches_compatible("001", "110"));
    BOOST_CHECK(TreeNode::branches_compatible("011", "110"));
    BOOST_CHECK(TreeNode::branches_compatible("011", "111"));
    BOOST_CHECK(TreeNode::branches_compatible("011", "100"));
    BOOST_CHECK(TreeNode::branches_compatible("0011", "0011"));
    BOOST_CHECK(TreeNode::branches_compatible("0011", "1100"));
    BOOST_CHECK(TreeNode::branches_compatible("0011", "0111"));
    BOOST_CHECK(!TreeNode::branches_compatible("0011", "0110"));
}

