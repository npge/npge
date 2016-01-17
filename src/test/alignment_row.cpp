/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"

BOOST_AUTO_TEST_CASE (MapAlignmentRow_main) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 9, 1);
    AlignmentRow* row = new MapAlignmentRow("TGGTCCGAGA");
    f1.set_row(row);
    BOOST_CHECK(row->length() == 10);
    BOOST_CHECK(row->map_to_alignment(3) == 3);
    BOOST_CHECK(row->map_to_fragment(3) == 3);
    BOOST_CHECK(row->nearest_in_fragment(3) == 3);
}

BOOST_AUTO_TEST_CASE (MapAlignmentRow_row_local) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    // row must not be defined here
    // row must not be deleted after fragment
    Fragment f1(s1, 0, 9, 1);
    MapAlignmentRow row("TGGTCCGAGA", &f1);
}

BOOST_AUTO_TEST_CASE (MapAlignmentRow_detach) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    MapAlignmentRow row("TGGTCCGAGA");
    Fragment f1(s1, 0, 9, 1);
    f1.set_row(&row);
    f1.detach_row();
}

BOOST_AUTO_TEST_CASE (MapAlignmentRow_2) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    // ----------------------------------------------------012345678901234567
    Fragment f1(s1, 0, 9, 1);
    AlignmentRow* row = new MapAlignmentRow("--TGGT---CCGAG-A--");
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
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 0, 9, 1);
    AlignmentRow* row = new CompactAlignmentRow("TGGTCCGAGA");
    f1.set_row(row);
    BOOST_CHECK(row->length() == 10);
    BOOST_CHECK(row->map_to_alignment(3) == 3);
    BOOST_CHECK(row->map_to_fragment(3) == 3);
    BOOST_CHECK(row->nearest_in_fragment(3) == 3);
}

BOOST_AUTO_TEST_CASE (CompactAlignmentRow_2) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    // ----------------------------------------------------012345678901234567
    Fragment f1(s1, 0, 9, 1);
    AlignmentRow* row = new CompactAlignmentRow("--TGGT---CCGAG-A--");
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
    using namespace npge;
    std::string seq, aln;
    for (int n = 0; n < NUMBER_OF_GROUPS; n++) {
        seq += std::string(n, 'A');
        aln += std::string(n, '-');
        aln += std::string(n, 'A');
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
    using namespace npge;
    SequencePtr seq = boost::make_shared<InMemorySequence>("GC");
    Fragment f1(seq, 0, seq->size());
    CompactAlignmentRow row;
    f1.set_row(&row);
    row.grow("G-----------------------------");
    row.grow("-----------------------------C");
}

BOOST_AUTO_TEST_CASE (AlignmentRow_434) {
    using namespace npge;
    CompactAlignmentRow row;
    row.grow("GCTGAAGCTGCCTGCATCGGTCGCTCGCGCGGTGGATTGACGACCAAGCTGCATGCTGTT");
    row.grow("GTCGATGCTATCGGCCTACCGCTGCGAATAAAGCCAACACCCGGCCATTATGGTGACTGT");
    row.grow("CCGCAAGCTTCAAGCCTTCTATCCGGCTTGAAGGGTGTGGGGCATGTCATTGCTGATGCA");
    row.grow("GCCTATGATGCCGATCACTTAAGGGCCTTCATTGCCAGCGATCTCAAGGCAACGGCTCAG");
    row.grow("ATCAAGGTCAATCCAACACGTTCCAGTGTCCCAACAATCGACTGGAGGCTGTACAAGGAA");
    row.grow("CGCCATCAGATTGAATGCTTTTTTAACAAGTTGAAACGCTATCGTCGTATTGCGCTGCGA");
    row.grow("TGCGAGAAAACATTGACCGCATTCATGGGTTTCGTCCATCTCGCATGCGCTATGATCTGG");
    row.grow("TTACGTTAAATGCAG");
    BOOST_CHECK(row.length() == 435);
}

BOOST_AUTO_TEST_CASE (AlignmentRow_map_to_alignment) {
    using namespace npge;
    std::vector<RowType> types;
    types.push_back(MAP_ROW);
    types.push_back(COMPACT_ROW);
    for (int length = 0; length < 150; length++) {
        BOOST_FOREACH (RowType type, types) {
            AlignmentRow* row = AlignmentRow::new_row(type);
            row->grow(std::string(length, 'A'));
            BOOST_REQUIRE(row->length() == length);
            for (int i = -150; i < 0; i++) {
                BOOST_CHECK(row->map_to_alignment(i) == -1);
            }
            for (int i = length; i < length + 150; i++) {
                BOOST_CHECK(row->map_to_alignment(i) == -1);
            }
            delete row;
        }
    }
}

BOOST_AUTO_TEST_CASE (AlignmentRow_InversedRow) {
    using namespace npge;
    SequencePtr s((new InMemorySequence("AATG")));
    boost::scoped_ptr<Fragment> f((new Fragment(s, 0, 3)));
    f->set_row(new CompactAlignmentRow("A-ATG"));
    BOOST_CHECK(f->row()->length() == 5);
    f->inverse();
    BOOST_CHECK(f->str() == "CAT-T");
    f->inverse();
    BOOST_CHECK(f->str() == "A-ATG");
    f->set_ori(-1);
    BOOST_CHECK(f->str() == "CAT-T");
}

