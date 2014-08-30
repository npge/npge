/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "class_name.hpp"

using namespace npge;

BOOST_AUTO_TEST_CASE (class_name_main) {
    // gcc
    BOOST_CHECK(class_name("N12npge6FilterE") == "Filter");
    BOOST_CHECK(class_name("N12npge9InE") == "In");
    BOOST_CHECK(class_name("N12npge17OverlapsResolver2E") ==
                "OverlapsResolver2");
    // MSVC
    BOOST_CHECK(class_name("class OverlapsResolver2") == "OverlapsResolver2");
    BOOST_CHECK(class_name("class Filter") == "Filter");
}

