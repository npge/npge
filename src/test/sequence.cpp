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
    BOOST_REQUIRE(seq.approximate_size() == 36);
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgg");
    BOOST_REQUIRE(s1->approximate_size() == 3);
    size_t fragments_number = 0;
    Fragment f(s1);
    s1->make_first_fragment(f, 2);
    while (s1->next_fragment(f)) {
        fragments_number += 1;
    }
    BOOST_REQUIRE(fragments_number == 4);
}

