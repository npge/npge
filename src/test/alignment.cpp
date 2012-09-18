/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Alignment.hpp"

BOOST_AUTO_TEST_CASE (Alignment_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    Alignment a;
    BOOST_REQUIRE(a.add_fragment(&f1, "tggtccgaga") == 0);
    BOOST_CHECK(a.length() == 10);
    BOOST_CHECK(a.size() == 1);
    BOOST_CHECK(a.index_of(&f1) == 0);
    BOOST_CHECK(a.fragment_at(0) == &f1);
    BOOST_CHECK(a.map_to_alignment(0, 3) == 3);
    BOOST_CHECK(a.map_to_fragment(0, 3) == 3);
    BOOST_CHECK(a.nearest_in_fragment(0, 3) == 3);
    a.remove_fragment(0);
    // TODO BOOST_CHECK(a.length() == 0);
    BOOST_CHECK(a.size() == 0);
}

