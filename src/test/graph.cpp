/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Graph.hpp"

using namespace npge;

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
    using namespace npge;
    Graph<int> b;
    b.push_back(std::make_pair(0, 1));
    BOOST_CHECK(!b.is_symmetric());
    b.add_for_symmetric();
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

static void test_connected_component2(SortedVector<int>& vertices,
                                      Graph<int>& edges) {
    BOOST_CHECK(vertices.size() == 10);
    BOOST_CHECK(edges.size() == 9);
}

BOOST_AUTO_TEST_CASE (graph_clique) {
    using namespace npge;
    Graph<int> b;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            b.push_back(std::make_pair(i, j));
        }
    }
    BOOST_CHECK(b.is_symmetric());
    b.connected_components(test_connected_component2);
    b.remove_extra_edges();
    b.connected_components(test_connected_component2);
}

