/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "BloomFilter.hpp"

BOOST_AUTO_TEST_CASE (BloomFilter_test) {
    bloomrepeats::BloomFilter filter(1e6, 0.01);
    BOOST_REQUIRE(filter.bits() == 9585058);
    BOOST_REQUIRE(filter.hashes() == 7);
    filter.add("atgc");
    filter.add("aaaa");
    filter.add("tatg", -1);
    BOOST_REQUIRE(filter.test("atgc"));
    BOOST_REQUIRE(filter.test("gcat", -1));
    BOOST_REQUIRE(filter.test("aaaa"));
    BOOST_REQUIRE(filter.test("tttt", -1));
    BOOST_REQUIRE(filter.test("tatg", -1));
    BOOST_REQUIRE(filter.test("cata"));
    BOOST_CHECK(!filter.test("gggg"));
    BOOST_CHECK(!filter.test("ttgc"));
}

BOOST_AUTO_TEST_CASE (BloomFilter_default_constructor) {
    bloomrepeats::BloomFilter filter;
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

BOOST_AUTO_TEST_CASE (BloomFilter_set_members) {
    bloomrepeats::BloomFilter filter;
    filter.set_members(2, 0.000001);
    filter.set_optimal_hashes(2);
    BOOST_REQUIRE(filter.bits() == 58);
    BOOST_REQUIRE(filter.hashes() == 20);
    filter.add("atgc");
    filter.add("aaaa");
    BOOST_REQUIRE(filter.test("atgc"));
    BOOST_REQUIRE(filter.test("aaaa"));
    BOOST_CHECK(!filter.test("gggg"));
    BOOST_CHECK(!filter.test("ttgc"));
}

BOOST_AUTO_TEST_CASE (BloomFilter_test_and_add) {
    bloomrepeats::BloomFilter filter(1e6, 0.01);
    BOOST_REQUIRE(filter.test_and_add("atgc") == false);
    BOOST_CHECK(filter.test_and_add("aaaa") == false);
    filter.add("ttaa");
    BOOST_REQUIRE(filter.test_and_add("atgc") == true);
    BOOST_REQUIRE(filter.test("atgc") == true);
    BOOST_REQUIRE(filter.test_and_add("ttaa") == true);
    BOOST_REQUIRE(filter.test("ttaa") == true);
}

