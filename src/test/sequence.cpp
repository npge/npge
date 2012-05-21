/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"

BOOST_AUTO_TEST_CASE (Sequence_main) {
    using namespace bloomrepeats;
    InMemorySequence seq("tggtccgagatgcgggcccgtaagcttacatacagg");
    BOOST_REQUIRE(seq.size() == 36);
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgg");
    BOOST_REQUIRE(s1->size() == 3);
    size_t fragments_number = 0;
    Fragment f(s1);
    s1->make_first_fragment(f, 2);
    while (s1->next_fragment(f)) {
        fragments_number += 1;
    }
    BOOST_CHECK(fragments_number == 4);
}

BOOST_AUTO_TEST_CASE (Sequence_first_ori) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgg");
    Fragment f(s1);
    s1->make_first_fragment(f, 2);
    while (s1->next_fragment(f)) {
        BOOST_CHECK(f.ori() == Sequence::FIRST_ORI);
        break;
    }
}

BOOST_AUTO_TEST_CASE (Sequence_filtering) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>(" ---ATG--caggacg..");
    BOOST_REQUIRE(s1->size() == 10);
    Fragment f(s1, 0, 9);
    BOOST_CHECK(f.str() == "atgcaggacg");
    BOOST_CHECK(f.at(-1) == 'g');
}

