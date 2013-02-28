/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"

BOOST_AUTO_TEST_CASE (MapAlignmentRow_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    AlignmentRow* row = new MapAlignmentRow("tggtccgaga");
    f1.set_row(row);
    BOOST_CHECK(row->length() == 10);
    BOOST_CHECK(row->map_to_alignment(3) == 3);
    BOOST_CHECK(row->map_to_fragment(3) == 3);
    BOOST_CHECK(row->nearest_in_fragment(3) == 3);
}

BOOST_AUTO_TEST_CASE (MapAlignmentRow_row_local) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    // row must not be defined here
    // row must not be deleted after fragment
    Fragment f1(s1, 0, 9, 1);
    MapAlignmentRow row("tggtccgaga", &f1);
}

BOOST_AUTO_TEST_CASE (MapAlignmentRow_detach) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    MapAlignmentRow row("tggtccgaga");
    Fragment f1(s1, 0, 9, 1);
    f1.set_row(&row);
    f1.detach_row();
}

BOOST_AUTO_TEST_CASE (MapAlignmentRow_2) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    // ----------------------------------------------------012345678901234567
    Fragment f1(s1, 0, 9, 1);
    AlignmentRow* row = new MapAlignmentRow("--tggt---ccgag-a--");
    // --------------------------------------012345678901234567
    f1.set_row(row);
    BOOST_CHECK(row->length() == 18);
    BOOST_CHECK(row->map_to_alignment(0) == 2);
    BOOST_CHECK(row->map_to_fragment(2) == 0);
    BOOST_CHECK(row->map_to_alignment(8) == 13);
    BOOST_CHECK(row->map_to_fragment(13) == 8);
    BOOST_CHECK(row->nearest_in_fragment(0) == 0);
    BOOST_CHECK(row->nearest_in_fragment(6) == 3);
}

BOOST_AUTO_TEST_CASE (CompactAlignmentRow_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 0, 9, 1);
    AlignmentRow* row = new CompactAlignmentRow("tggtccgaga");
    f1.set_row(row);
    BOOST_CHECK(row->length() == 10);
    BOOST_CHECK(row->map_to_alignment(3) == 3);
    BOOST_CHECK(row->map_to_fragment(3) == 3);
    BOOST_CHECK(row->nearest_in_fragment(3) == 3);
}

BOOST_AUTO_TEST_CASE (CompactAlignmentRow_2) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    // ----------------------------------------------------012345678901234567
    Fragment f1(s1, 0, 9, 1);
    AlignmentRow* row = new CompactAlignmentRow("--tggt---ccgag-a--");
    // ------------------------------------------012345678901234567
    f1.set_row(row);
    BOOST_CHECK(row->length() == 18);
    BOOST_CHECK(row->map_to_alignment(0) == 2);
    BOOST_CHECK(row->map_to_fragment(2) == 0);
    BOOST_CHECK(row->map_to_alignment(8) == 13);
    BOOST_CHECK(row->map_to_fragment(13) == 8);
    BOOST_CHECK(row->nearest_in_fragment(0) == 0);
    BOOST_CHECK(row->nearest_in_fragment(6) == 3);
}

BOOST_AUTO_TEST_CASE (CompactAlignmentRow_3) {
    const int NUMBER_OF_GROUPS = 100;
    using namespace bloomrepeats;
    std::string seq, aln;
    for (int n = 0; n < NUMBER_OF_GROUPS; n++) {
        seq += std::string(n, 'a');
        aln += std::string(n, '-');
        aln += std::string(n, 'a');
    }
    SequencePtr s1 = boost::make_shared<InMemorySequence>(seq);
    Fragment f1(s1, 0, seq.length());
    AlignmentRow* row = new CompactAlignmentRow(aln);
    f1.set_row(row);
    BOOST_CHECK(row->length() == aln.length());
    std::stringstream row_str;
    f1.print_contents(row_str, '-', /* line */ 0);
    BOOST_CHECK(row_str.str() == aln);
    for (int align_pos = 0; align_pos < aln.length(); align_pos++) {
        int f_pos = row->map_to_fragment(align_pos);
        if (f_pos != -1) {
            BOOST_CHECK(row->map_to_alignment(f_pos) == align_pos);
        }
    }
}

BOOST_AUTO_TEST_CASE (AlignmentRow_letter_and_gaps_in_begin) {
    using namespace bloomrepeats;
    SequencePtr seq = boost::make_shared<InMemorySequence>("gc");
    Fragment f1(seq, 0, seq->size());
    CompactAlignmentRow row;
    f1.set_row(&row);
    row.grow("g-----------------------------");
    row.grow("-----------------------------c");
}

