/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"

BOOST_AUTO_TEST_CASE (Fragment_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    BOOST_REQUIRE(f1.length() == 10);
    BOOST_REQUIRE(f1.str() == "tggtccgaga");
    BOOST_REQUIRE(*f1.begin() == 't');
    Fragment f2(s1, 0, 9, -1);
    BOOST_REQUIRE(f2.length() == 10);
    BOOST_REQUIRE(f2.str() == "tctcggacca");
    BOOST_REQUIRE(*f2.begin() == 'a');
    BOOST_REQUIRE(f2.begin() - f2.end() == 10);
}

BOOST_AUTO_TEST_CASE (Fragment_expand) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagaTgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    BOOST_REQUIRE(f1.length() == 10);
    BOOST_REQUIRE(f1.expand() == 't');
    BOOST_REQUIRE(f1.length() == 11);
    f1.inverse();
    BOOST_REQUIRE(!f1.expand());
    BOOST_REQUIRE(f1.length() == 11);
}

