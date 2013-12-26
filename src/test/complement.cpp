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
    BOOST_CHECK(complement('A') == 'T');
    BOOST_CHECK(complement('T') == 'A');
    BOOST_CHECK(complement('G') == 'C');
    BOOST_CHECK(complement('C') == 'G');
}

BOOST_AUTO_TEST_CASE (complement_string) {
    using namespace bloomrepeats;
    std::string data;
    data = "AAA";
    complement(data);
    BOOST_CHECK(data == "TTT");
    //
    data = "ATGC";
    complement(data);
    BOOST_CHECK(data == "GCAT");
    //
    data = "ACGT";
    complement(data);
    BOOST_CHECK(data == "ACGT");
    //
    data = "A-CG-T";
    complement(data);
    BOOST_CHECK(data == "A-CG-T");
    //
    data = "TT~~A";
    complement(data);
    BOOST_CHECK(data == "T~~AA");
}

