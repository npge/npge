/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"

BOOST_AUTO_TEST_CASE (MapAlignmentRow_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    MapAlignmentRow row(&f1, "tggtccgaga");
    BOOST_CHECK(row.length() == 10);
    BOOST_CHECK(row.map_to_alignment(3) == 3);
    BOOST_CHECK(row.map_to_fragment(3) == 3);
    BOOST_CHECK(row.nearest_in_fragment(3) == 3);
}

BOOST_AUTO_TEST_CASE (MapAlignmentRow_2) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    // ----------------------------------------------------012345678901234567
    Fragment f1(s1, 0, 9, 1);
    MapAlignmentRow row(&f1, "--tggt---ccgag-a--");
    // -----------------------012345678901234567
    BOOST_CHECK(row.length() == 18);
    BOOST_CHECK(row.map_to_alignment(0) == 2);
    BOOST_CHECK(row.map_to_fragment(2) == 0);
    BOOST_CHECK(row.map_to_alignment(8) == 13);
    BOOST_CHECK(row.map_to_fragment(13) == 8);
    BOOST_CHECK(row.nearest_in_fragment(0) == 0);
    BOOST_CHECK(row.nearest_in_fragment(6) == 3);
}

BOOST_AUTO_TEST_CASE (CompactAlignmentRow_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    CompactAlignmentRow row(&f1, "tggtccgaga");
    BOOST_CHECK(row.length() == 10);
    BOOST_CHECK(row.map_to_alignment(3) == 3);
    BOOST_CHECK(row.map_to_fragment(3) == 3);
    BOOST_CHECK(row.nearest_in_fragment(3) == 3);
}

BOOST_AUTO_TEST_CASE (CompactAlignmentRow_2) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    // ----------------------------------------------------012345678901234567
    Fragment f1(s1, 0, 9, 1);
    CompactAlignmentRow row(&f1, "--tggt---ccgag-a--");
    // ---------------------------012345678901234567
    BOOST_CHECK(row.length() == 18);
    BOOST_CHECK(row.map_to_alignment(0) == 2);
    BOOST_CHECK(row.map_to_fragment(2) == 0);
    BOOST_CHECK(row.map_to_alignment(8) == 13);
    BOOST_CHECK(row.map_to_fragment(13) == 8);
    BOOST_CHECK(row.nearest_in_fragment(0) == 0);
    BOOST_CHECK(row.nearest_in_fragment(6) == 3);
}

