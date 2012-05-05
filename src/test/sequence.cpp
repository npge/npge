/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"

BOOST_AUTO_TEST_CASE (Sequence_main) {
    bloomrepeats::InMemorySequence seq("tggtccgagatgcgggcccgtaagcttacatacagg");
    BOOST_REQUIRE(seq.approximate_size() == 36);
    size_t length = 5;
    BOOST_REQUIRE(*seq.get(7, length) == 'a');
    BOOST_REQUIRE(length == 5);
    BOOST_REQUIRE(*seq.get(34, length) == 'g');
    BOOST_REQUIRE(length == 2);
}

