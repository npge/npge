/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "class_name.hpp"

using namespace bloomrepeats;

BOOST_AUTO_TEST_CASE (class_name_main) {
    // gcc
    BOOST_CHECK(class_name("N12bloomrepeats6FilterE") == "Filter");
    BOOST_CHECK(class_name("N12bloomrepeats9AddBlocksE") == "AddBlocks");
    BOOST_CHECK(class_name("N12bloomrepeats17OverlapsResolver2E") ==
            "OverlapsResolver2");
    // MSVC
    BOOST_CHECK(class_name("class OverlapsResolver2") == "OverlapsResolver2");
    BOOST_CHECK(class_name("class Filter") == "Filter");
}

