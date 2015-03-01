/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "BloomFilter.hpp"

BOOST_AUTO_TEST_CASE (BloomFilter_test) {
    npge::BloomFilter filter(1e6, 0.01);
    BOOST_CHECK(filter.bits() == 9585059);
    BOOST_CHECK(filter.hashes() == 7);
    filter.add("ATGC");
    filter.add("AAAA");
    filter.add("TATG", -1);
    BOOST_CHECK(filter.test("ATGC"));
    BOOST_CHECK(filter.test("GCAT", -1));
    BOOST_CHECK(filter.test("AAAA"));
    BOOST_CHECK(filter.test("TTTT", -1));
    BOOST_CHECK(filter.test("TATG", -1));
    BOOST_CHECK(filter.test("CATA"));
    BOOST_WARN(!filter.test("GGGG"));
    BOOST_WARN(!filter.test("TTGC"));
}

BOOST_AUTO_TEST_CASE (BloomFilter_optimal) {
    using namespace npge;
    BOOST_CHECK(BloomFilter::optimal_bits(0, 0.1) < 100);
    BOOST_CHECK(BloomFilter::optimal_bits(0, 0.1) > 0);
    BOOST_CHECK(BloomFilter::optimal_bits(1, 0.1) < 100);
    BOOST_CHECK(BloomFilter::optimal_bits(1, 0.1) > 0);
    BOOST_CHECK(BloomFilter::optimal_hashes(0, 1) < 100);
    BOOST_CHECK(BloomFilter::optimal_hashes(0, 1) > 0);
    BOOST_CHECK(BloomFilter::optimal_hashes(1, 1) < 100);
    BOOST_CHECK(BloomFilter::optimal_hashes(1, 1) > 0);
}

BOOST_AUTO_TEST_CASE (BloomFilter_default_constructor) {
    npge::BloomFilter filter;
    filter.set_bits(50);
    filter.set_hashes(5);
    BOOST_REQUIRE(filter.bits() == 50);
    BOOST_REQUIRE(filter.hashes() == 5);
    filter.add("ATGC");
    filter.add("AAAA");
    BOOST_CHECK(filter.test("ATGC"));
    BOOST_CHECK(filter.test("AAAA"));
    BOOST_WARN(!filter.test("GGGG"));
    BOOST_WARN(!filter.test("TTGC"));
}

BOOST_AUTO_TEST_CASE (BloomFilter_set_members) {
    npge::BloomFilter filter;
    filter.set_members(2, 0.000001);
    filter.set_optimal_hashes(2);
    BOOST_CHECK(filter.bits() == 59);
    BOOST_CHECK(filter.hashes() == 20);
    filter.add("ATGC");
    filter.add("AAAA");
    BOOST_CHECK(filter.test("ATGC"));
    BOOST_CHECK(filter.test("AAAA"));
    BOOST_WARN(!filter.test("GGGG"));
    BOOST_WARN(!filter.test("TTGC"));
}

BOOST_AUTO_TEST_CASE (BloomFilter_test_and_add) {
    npge::BloomFilter filter(1e6, 0.01);
    BOOST_CHECK(filter.test_and_add("ATGC") == false);
    BOOST_WARN(filter.test_and_add("AAAA") == false);
    filter.add("TTAA");
    BOOST_CHECK(filter.test_and_add("ATGC") == true);
    BOOST_CHECK(filter.test("ATGC") == true);
    BOOST_CHECK(filter.test_and_add("TTAA") == true);
    BOOST_CHECK(filter.test("TTAA") == true);
}

