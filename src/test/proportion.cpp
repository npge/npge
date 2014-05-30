/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "proportion.hpp"

BOOST_AUTO_TEST_CASE (proportion_main) {
    using namespace npge;
    BOOST_CHECK(proportion(0, 5, 5) == 0);
    BOOST_CHECK(proportion(1, 5, 5) == 1);
    BOOST_CHECK(proportion(2, 5, 5) == 2);
    BOOST_CHECK(proportion(3, 5, 5) == 3);
    BOOST_CHECK(proportion(4, 5, 5) == 4);
    BOOST_CHECK(proportion(5, 5, 5) == 5);
    BOOST_CHECK(proportion(6, 5, 5) == 6);
    BOOST_CHECK(proportion(1, 6, 10) == 1);
    BOOST_CHECK(proportion(5, 10, 20) == 10);
    BOOST_CHECK(proportion(0, 10, 20) == 0);
    BOOST_CHECK(proportion(5, 10, 0) == 0);
    BOOST_CHECK(proportion(0, 0, 0) == 0);
    BOOST_CHECK(proportion(0, 10, 10) == 0);
    BOOST_CHECK(proportion(0, 0, 10) == 0);
    BOOST_CHECK(proportion(0, 10, 0) == 0);
    BOOST_CHECK(proportion(10, 10, 0) == 0);
    BOOST_CHECK(proportion(5e7, 1e8, 2e8) == 1e8);
    BOOST_CHECK(proportion(1, 2, 2e8) == 1e8);
    BOOST_CHECK(proportion(1, 1e8, 2e8) == 2);
    BOOST_CHECK(proportion(0, 1e8, 2e8) == 0);
    BOOST_CHECK(proportion(10, 1, 20) == 200);
    BOOST_CHECK(proportion(10, 20, 8) == 4);
    BOOST_CHECK(proportion(1e8, 2e8, 5e7) == 25e6);
    BOOST_CHECK(proportion(100000, 100001, 64) == 63);
}

