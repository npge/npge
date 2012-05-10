/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "complement.hpp"

BOOST_AUTO_TEST_CASE (complement_char) {
    using namespace bloomrepeats;
    BOOST_REQUIRE(complement('a') == 't');
    BOOST_REQUIRE(complement('t') == 'a');
    BOOST_REQUIRE(complement('g') == 'c');
    BOOST_REQUIRE(complement('c') == 'g');
}

BOOST_AUTO_TEST_CASE (complement_string) {
    using namespace bloomrepeats;
    std::string data;
    data = "aaa";
    complement(data);
    BOOST_REQUIRE(data == "ttt");
    //
    data = "atgc";
    complement(data);
    BOOST_REQUIRE(data == "gcat");
    //
    data = "acgt";
    complement(data);
    BOOST_REQUIRE(data == "acgt");
    //
    data = "a-cg-t";
    complement(data);
    BOOST_REQUIRE(data == "a-cg-t");
    //
    data = "tt~~a";
    complement(data);
    BOOST_REQUIRE(data == "t~~aa");
}

