/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Graph.hpp"

using namespace bloomrepeats;

static void test_connected_component(SortedVector<int>& vertices,
                                     Graph<int>& edges) {
    BOOST_CHECK(vertices.size() == 4);
    BOOST_CHECK(edges.size() == 3);
    vertices.sort_unique();
    edges.sort_unique();
    BOOST_CHECK(vertices.size() == 4);
    BOOST_CHECK(edges.size() == 3);
}

BOOST_AUTO_TEST_CASE (graph_main) {
    using namespace bloomrepeats;
    Graph<int> b;
    b.push_back(std::make_pair(0, 1));
    BOOST_CHECK(!b.is_symmetric());
    b.push_back(std::make_pair(1, 0));
    BOOST_CHECK(b.is_symmetric());
    b.push_back(std::make_pair(0, 2));
    b.push_back(std::make_pair(2, 0));
    b.push_back(std::make_pair(1, 3));
    b.push_back(std::make_pair(3, 1));
    b.sort_unique();
    BOOST_CHECK(b.is_symmetric());
    SortedVector<int> v;
    b.connected_with(v, 0);
    BOOST_CHECK(v.size() == 2 && v[0] == 1 && v[1] == 2);
    v.clear();
    b.connected_with(v, 1);
    BOOST_CHECK(v.size() == 2 && v[0] == 0 && v[1] == 3);
    v.clear();
    b.connected_with(v, 2);
    BOOST_CHECK(v.size() == 1 && v[0] == 0);
    v.clear();
    b.connected_with(v, 3);
    BOOST_CHECK(v.size() == 1 && v[0] == 1);
    //
    v.clear();
    b.all_vertices(v);
    BOOST_CHECK(v.is_sorted_unique());
    BOOST_CHECK(v.size() == 4);
    b.connected_components(test_connected_component);
}

