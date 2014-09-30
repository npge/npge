/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "boundaries.hpp"

BOOST_AUTO_TEST_CASE (boundaries_avg_element) {
    using namespace npge;
    Boundaries b;
    b.push_back(0);
    BOOST_CHECK(avg_element(b) == 0);
    b.push_back(2);
    BOOST_CHECK(avg_element(b) == 1);
    Boundaries a;
    a.push_back(1e9);
    a.push_back(1e9);
    a.push_back(1e9);
    a.push_back(1e9);
    BOOST_CHECK(avg_element(a) == 1e9);
}

BOOST_AUTO_TEST_CASE (boundaries_nearest_element) {
    using namespace npge;
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
    b.push_back(1e9);
    b.push_back(2e9);
    BOOST_CHECK(nearest_element(b, 100) == 110);
    b.push_back(2e9 + 1);
    BOOST_CHECK(nearest_element(b, 2e9) == 2e9);
    BOOST_CHECK(nearest_element(b, 2e9 + 100) == 2e9 + 1);
    Boundaries a;
    a.push_back(5);
    a.push_back(15);
    BOOST_CHECK(nearest_element(a, 10) == 15); // to largest
}

BOOST_AUTO_TEST_CASE (boundaries_bound) {
    using namespace npge;
    Boundaries b;
    b.push_back(0);
    b.push_back(1);
    b.push_back(110);
    b.push_back(1e9);
    b.push_back(2e9);
    b.push_back(2e9 + 1);
    BOOST_CHECK(b.lower_bound(0) == b.begin());
    BOOST_CHECK(*b.lower_bound(1) == 1);
    BOOST_CHECK(*b.lower_bound(2) == 110);
    BOOST_CHECK(*b.lower_bound(2e9) == 2e9);
    BOOST_CHECK(b.lower_bound(2e9 + 10) == b.end());
    BOOST_CHECK(*b.upper_bound(0) == 1);
    BOOST_CHECK(*b.upper_bound(110) == 1e9);
    BOOST_CHECK(*b.upper_bound(2e9) == 2e9 + 1);
    BOOST_CHECK(b.has_elem(0));
    BOOST_CHECK(b.has_elem(1));
    BOOST_CHECK(!b.has_elem(2));
    BOOST_CHECK(b.has_elem(110));
    BOOST_CHECK(b.has_elem(2e9 + 1));
    BOOST_CHECK(!b.has_elem(2e9 - 1));
}

