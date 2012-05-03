/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "BloomFilter.hpp"

BOOST_AUTO_TEST_CASE (BloomFilter_test) {
    br::BloomFilter filter(100, 2);
    BOOST_REQUIRE(filter.bits() == 100);
    BOOST_REQUIRE(filter.hashes() == 2);
    filter.add("atgc");
    filter.add("aaaa");
    BOOST_REQUIRE(filter.test("atgc"));
    BOOST_REQUIRE(filter.test("aaaa"));
    BOOST_CHECK(!filter.test("gggg"));
    BOOST_CHECK(!filter.test("ttgc"));
}

BOOST_AUTO_TEST_CASE (BloomFilter_default_constructor) {
    br::BloomFilter filter;
    filter.set_bits(50);
    filter.set_hashes(5);
    BOOST_REQUIRE(filter.bits() == 50);
    BOOST_REQUIRE(filter.hashes() == 5);
    filter.add("atgc");
    filter.add("aaaa");
    BOOST_REQUIRE(filter.test("atgc"));
    BOOST_REQUIRE(filter.test("aaaa"));
    BOOST_CHECK(!filter.test("gggg"));
    BOOST_CHECK(!filter.test("ttgc"));
}

