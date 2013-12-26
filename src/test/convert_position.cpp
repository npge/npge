/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "convert_position.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Block.hpp"

BOOST_AUTO_TEST_CASE (convert_pos_alignment) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Fragment* f1 = new Fragment(s1, 2, 6, 1); // long
    f1->set_row(new CompactAlignmentRow("GT-CCG"));
    Fragment* f2 = new Fragment(s1, 9, 10, -1); // short
    f2->set_row(new CompactAlignmentRow("--G-T-"));
    Block block;
    block.insert(f1);
    block.insert(f2);
    int block_length = block.alignment_length();
    BOOST_REQUIRE(block_length == 6);
    BOOST_CHECK(block_pos(f1, 0, block_length) == 0);
    BOOST_CHECK(fragment_pos(f1, 0, block_length) == 0);
    BOOST_CHECK(block_pos(f1, 1, block_length) == 1);
    BOOST_CHECK(fragment_pos(f1, 1, block_length) == 1);
    BOOST_CHECK(fragment_pos(f1, 2, block_length) == 1 ||
                fragment_pos(f1, 2, block_length) == 2);
    BOOST_CHECK(block_pos(f1, 2, block_length) == 3);
    BOOST_CHECK(fragment_pos(f1, 3, block_length) == 2);
    BOOST_CHECK(block_pos(f1, 3, block_length) == 4);
    BOOST_CHECK(fragment_pos(f1, 4, block_length) == 3);
    BOOST_CHECK(block_pos(f1, 4, block_length) == 5);
    BOOST_CHECK(fragment_pos(f1, 5, block_length) == 4);
    //
    BOOST_CHECK(fragment_pos(f2, 0, block_length) == 0);
    BOOST_CHECK(fragment_pos(f2, 1, block_length) == 0);
    BOOST_CHECK(block_pos(f2, 0, block_length) == 2);
    BOOST_CHECK(fragment_pos(f2, 2, block_length) == 0);
    BOOST_CHECK(fragment_pos(f2, 3, block_length) == 0 ||
                fragment_pos(f2, 3, block_length) == 1);
    BOOST_CHECK(block_pos(f2, 1, block_length) == 4);
    BOOST_CHECK(fragment_pos(f2, 4, block_length) == 1);
    BOOST_CHECK(fragment_pos(f2, 5, block_length) == 1);
}

BOOST_AUTO_TEST_CASE (convert_pos_noalignment) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Fragment* f1 = new Fragment(s1, 2, 6, 1); // long
    Fragment* f2 = new Fragment(s1, 9, 10, -1); // short
    Block block;
    block.insert(f1);
    block.insert(f2);
    int block_length = block.alignment_length();
    BOOST_REQUIRE(block_length == 5);
    BOOST_CHECK(block_pos(f1, 0, block_length) == 0);
    BOOST_CHECK(block_pos(f1, 1, block_length) == 1);
    BOOST_CHECK(block_pos(f1, 2, block_length) == 2);
    BOOST_CHECK(block_pos(f1, 3, block_length) == 3);
    BOOST_CHECK(block_pos(f1, 4, block_length) == 4);
    BOOST_CHECK(block_pos(f2, 0, block_length) == 0);
    BOOST_CHECK(block_pos(f2, 1, block_length) == 2);
}

BOOST_AUTO_TEST_CASE (convert_seq) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Fragment f1(s1, 2, 6);
    BOOST_CHECK(frag_to_seq(&f1, 0) == 2);
    BOOST_CHECK(seq_to_frag(&f1, 2) == 0);
    BOOST_CHECK(frag_to_seq(&f1, 1) == 3);
    BOOST_CHECK(seq_to_frag(&f1, 3) == 1);
    BOOST_CHECK(frag_to_seq(&f1, -1) == 1);
    BOOST_CHECK(seq_to_frag(&f1, 1) == -1);
    Fragment f2(s1, 2, 6, -1);
    BOOST_CHECK(frag_to_seq(&f2, 0) == 6);
    BOOST_CHECK(seq_to_frag(&f2, 6) == 0);
    BOOST_CHECK(frag_to_seq(&f2, 1) == 5);
    BOOST_CHECK(seq_to_frag(&f2, 5) == 1);
    BOOST_CHECK(frag_to_seq(&f2, -1) == 7);
    BOOST_CHECK(seq_to_frag(&f2, 7) == -1);
}

