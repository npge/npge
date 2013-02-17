/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "boundaries.hpp"

BOOST_AUTO_TEST_CASE (boundaries_avg_element) {
    using namespace bloomrepeats;
    Boundaries b;
    b.push_back(0);
    BOOST_CHECK(avg_element(b) == 0);
    b.push_back(2);
    BOOST_CHECK(avg_element(b) == 1);
    Boundaries a;
    a.push_back(3e9);
    a.push_back(3e9);
    a.push_back(3e9);
    a.push_back(3e9);
    BOOST_CHECK(avg_element(a) == 3e9);
}

BOOST_AUTO_TEST_CASE (boundaries_nearest_element) {
    using namespace bloomrepeats;
    Boundaries b;
    b.push_back(0);
    BOOST_CHECK(nearest_element(b, 0) == 0);
    BOOST_CHECK(nearest_element(b, 1) == 0);
    b.push_back(1);
    BOOST_CHECK(nearest_element(b, 0) == 0);
    BOOST_CHECK(nearest_element(b, 1) == 1);
    BOOST_CHECK(nearest_element(b, 100) == 1);
    b.push_back(110);
    BOOST_CHECK(nearest_element(b, 100) == 110);
    b.push_back(2e9);
    b.push_back(3e9);
    BOOST_CHECK(nearest_element(b, 100) == 110);
    b.push_back(3e9 + 1);
    BOOST_CHECK(nearest_element(b, 3e9) == 3e9);
    BOOST_CHECK(nearest_element(b, 4e9) == 3e9 + 1);
    BOOST_CHECK(nearest_element(b, 21e8) == 2e9);
    BOOST_CHECK(nearest_element(b, 29e8) == 3e9);
    BOOST_CHECK(nearest_element(b, 25e8) == 3e9);
    Boundaries a;
    a.push_back(5);
    a.push_back(15);
    BOOST_CHECK(nearest_element(a, 10) == 15); // to largest
}

BOOST_AUTO_TEST_CASE (boundaries_select_boundaries) {
    using namespace bloomrepeats;
    Boundaries b;
    b.push_back(100);
    b.push_back(100);
    b.push_back(101);
    b.push_back(3e9);
    b.push_back(3e9 + 10);
    b.push_back(3e9 + 100);
    b.push_back(0);
    b.push_back(1);
    b.push_back(5);
    select_boundaries(b, 20, 3e9 + 110);
    BOOST_CHECK(b.size() == 4);
    BOOST_CHECK(b[0] == 0);
    BOOST_CHECK(b[1] == 100);
    BOOST_CHECK(b[2] == 3e9 + 5 && b[3] == 3e9 + 110);
}

BOOST_AUTO_TEST_CASE (boundaries_bound) {
    using namespace bloomrepeats;
    Boundaries b;
    b.push_back(0);
    b.push_back(1);
    b.push_back(110);
    b.push_back(2e9);
    b.push_back(3e9);
    b.push_back(3e9 + 1);
    BOOST_CHECK(b.lower_bound(0) == b.begin());
    BOOST_CHECK(*b.lower_bound(1) == 1);
    BOOST_CHECK(*b.lower_bound(2) == 110);
    BOOST_CHECK(*b.lower_bound(3e9) == 3e9);
    BOOST_CHECK(b.lower_bound(3e9 + 10) == b.end());
    BOOST_CHECK(*b.upper_bound(0) == 1);
    BOOST_CHECK(*b.upper_bound(110) == 2e9);
    BOOST_CHECK(*b.upper_bound(3e9) == 3e9 + 1);
    BOOST_CHECK(b.has_elem(0));
    BOOST_CHECK(b.has_elem(1));
    BOOST_CHECK(!b.has_elem(2));
    BOOST_CHECK(b.has_elem(110));
    BOOST_CHECK(b.has_elem(3e9 + 1));
    BOOST_CHECK(!b.has_elem(3e9 - 1));
}

