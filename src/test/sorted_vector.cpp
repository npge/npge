/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/test/unit_test.hpp>

#include "SortedVector.hpp"

BOOST_AUTO_TEST_CASE (sorted_vector_main) {
    using namespace bloomrepeats;
    SortedVector<int> b;
    b.push_back(0);
    b.push_back(5);
    b.push_back(1);
    BOOST_CHECK(!b.is_sorted_unique());
    b.sort_unique();
    BOOST_CHECK(b.is_sorted_unique());
    BOOST_CHECK(b[0] == 0 && b[1] == 1 && b[2] == 5);
    b.push_back(5);
    BOOST_CHECK(!b.is_sorted_unique());
    b.sort_unique();
    BOOST_CHECK(b.is_sorted_unique());
    BOOST_CHECK(b.size() == 3 && b[0] == 0 && b[1] == 1 && b[2] == 5);
    BOOST_CHECK(b.has_elem(0));
    BOOST_CHECK(b.has_elem(1));
    BOOST_CHECK(b.has_elem(5));
    BOOST_CHECK(!b.has_elem(2));
    BOOST_CHECK(!b.has_elem(6));
    SortedVector<int> a;
    a.extend(b);
    a.extend(b);
    BOOST_CHECK(a.size() == 2 * b.size());
    std::stringstream ss;
    ss << b;
    BOOST_CHECK(ss.str() == "0\n1\n5\n");
}

