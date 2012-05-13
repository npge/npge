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
    BOOST_REQUIRE(f1.substr(1, 1) == "g");
    BOOST_REQUIRE(f1.substr(1, 2) == "gg");
    BOOST_REQUIRE(f1.substr(1, -1) == "ggtccgaga");
    BOOST_REQUIRE(f1.substr(-2, -1) == "ga");
    BOOST_REQUIRE(f2.substr(1, 1) == "c");
    BOOST_REQUIRE(f2.substr(1, 2) == "ct");
    BOOST_REQUIRE(f2.substr(1, -1) == "ctcggacca");
    BOOST_REQUIRE(f2.substr(-2, -1) == "ca");
}

BOOST_AUTO_TEST_CASE (Fragment_expand) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagaTgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    BOOST_REQUIRE(f1.length() == 10);
    f1.shift_end();
    BOOST_REQUIRE(f1.valid());
    BOOST_REQUIRE(f1.length() == 11);
    f1.inverse();
    f1.shift_end();
    BOOST_REQUIRE(!f1.valid());
    f1.shift_end(-1);
    BOOST_REQUIRE(f1.length() == 11);
    BOOST_REQUIRE(f1.valid());
    f1.shift_end(-1);
    BOOST_REQUIRE(f1.length() == 10);
    BOOST_REQUIRE(f1.min_pos() == 1);
}

BOOST_AUTO_TEST_CASE (Fragment_equal) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    Fragment f2(s1, 0, 9, 1);
    Fragment f3(s1, 0, 9, -1);
    BOOST_REQUIRE(f1 == f1);
    BOOST_REQUIRE(f1 == f2);
    BOOST_REQUIRE(f1 != f3);
}

BOOST_AUTO_TEST_CASE (Fragment_at) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment f1(s1, 5, 10, 1);
    Fragment f2(s1, 5, 10, -1);
    BOOST_REQUIRE(f1.at(0) == 'c');
    BOOST_REQUIRE(f1.at(1) == 'g');
    BOOST_REQUIRE(f1.at(-1) == 't');
    BOOST_REQUIRE(f1.at(-2) == 'a');
    BOOST_REQUIRE(f2.at(0) == 'a');
    BOOST_REQUIRE(f2.at(1) == 't');
    BOOST_REQUIRE(f2.at(-1) == 'g');
    BOOST_REQUIRE(f2.at(-2) == 'c');
}

BOOST_AUTO_TEST_CASE (Fragment_next) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 1, 2, 1);
    FragmentPtr f2 = boost::make_shared<Fragment>(s1, 5, 6, -1);
    FragmentPtr f3 = boost::make_shared<Fragment>(s1, 7, 8, 1);
    Fragment::connect(f1, f2);
    Fragment::connect(f2, f3);
    Fragment::connect(f3, f1);
    BOOST_REQUIRE(f1->next() == f2);
    BOOST_REQUIRE(f2->next() == f3);
    BOOST_REQUIRE(f3->next() == f1);
    BOOST_REQUIRE(f1->prev() == f3);
    BOOST_REQUIRE(f2->prev() == f1);
    BOOST_REQUIRE(f3->prev() == f2);
    f2.reset();
    BOOST_REQUIRE(!f1->next());
    BOOST_REQUIRE(!f3->prev());
    BOOST_REQUIRE(f1->prev() == f3);
    BOOST_REQUIRE(f3->next() == f1);
    f1->disconnect();
    BOOST_REQUIRE(!f1->prev());
    BOOST_REQUIRE(!f1->next());
    BOOST_REQUIRE(!f3->prev());
    BOOST_REQUIRE(!f3->next());
}

