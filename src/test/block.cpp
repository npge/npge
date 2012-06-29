/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <vector>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "PairAligner.hpp"

BOOST_AUTO_TEST_CASE (Block_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    FragmentPtr f1 = Fragment::create_new(s1, 2, 6, 1);
    FragmentPtr f2 = Fragment::create_new(s1, 9, 13, -1);
    BOOST_REQUIRE(!f1->block());
    BlockPtr block = Block::create_new();
    BOOST_REQUIRE(!block->front());
    block->insert(f1);
    block->insert(f2);
    BOOST_REQUIRE(block->size() == 2);
    BOOST_REQUIRE(block->has(f1));
    BOOST_REQUIRE(block->has(f2));
    BOOST_CHECK(block->front() == f1 || block->front() == f2);
    BOOST_REQUIRE(f1->block() == block);
    BOOST_REQUIRE(f2->block() == block);
    BOOST_FOREACH (FragmentPtr fragment, *block) {
        BOOST_CHECK(fragment == f1 || fragment == f2);
    }
    block->erase(f1);
    BOOST_CHECK(block->size() == 1);
    BOOST_CHECK(f2->block());
    BOOST_CHECK(!block->has(f1));
    BOOST_CHECK(block->has(f2));
    BOOST_CHECK(block->front() == f2);
    block->clear();
    BOOST_CHECK(block->size() == 0);
    BOOST_CHECK(block->empty());
}

BOOST_AUTO_TEST_CASE (Block_identity) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtTcgagcAgTcggcc");
    BlockPtr b1 = Block::create_new();
    b1->insert(Fragment::create_new(s1, 0, s1->size() - 1));
    b1->insert(Fragment::create_new(s2, 0, s1->size() - 1));
    BOOST_CHECK(std::abs(b1->identity() - 15. / 18.) < 0.01);
    BlockPtr b2 = Block::create_new();
    b2->insert(Fragment::create_new(s1, 0, 2));
    b2->insert(Fragment::create_new(s2, 0, 5));
    BOOST_CHECK(std::abs(b2->identity() - 0.5) < 0.01);
    BlockPtr b3 = Block::create_new();
    b3->insert(Fragment::create_new(s1, 0, 2));
    b3->insert(Fragment::create_new(s2, 0, 2));
    BOOST_CHECK(std::abs(b3->identity() - 1) < 0.01);
}

BOOST_AUTO_TEST_CASE (Block_match) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    b1->insert(Fragment::create_new(s1, 1, 2));
    b1->insert(Fragment::create_new(s2, 1, 2));
    b2->insert(Fragment::create_new(s1, 5, 6));
    b2->insert(Fragment::create_new(s2, 5, 6));
    BOOST_CHECK(Block::match(b1, b2) == 1);
    b1->insert(Fragment::create_new(s1, 8, 9));
    BOOST_CHECK(!Block::match(b1, b2));
    b2->insert(Fragment::create_new(s1, 10, 11));
    BOOST_CHECK(Block::match(b1, b2) == 1);
    b1->insert(Fragment::create_new(s1, 12, 13));
    b2->insert(Fragment::create_new(s1, 14, 15, -1));
    BOOST_CHECK(!Block::match(b1, b2));
}

BOOST_AUTO_TEST_CASE (Block_match_1) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    b1->insert(Fragment::create_new(s1, 1, 2));
    b1->insert(Fragment::create_new(s2, 1, 2));
    b2->insert(Fragment::create_new(s1, 5, 6, -1));
    b2->insert(Fragment::create_new(s2, 5, 6, -1));
    BOOST_CHECK(Block::match(b1, b2) == -1);
}

BOOST_AUTO_TEST_CASE (Block_join) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    FragmentPtr f21 = Fragment::create_new(s1, 3, 4);
    FragmentPtr f22 = Fragment::create_new(s2, 3, 4);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BOOST_CHECK(Block::can_join(b1, b2) == 1);
    BOOST_CHECK(Block::can_join(b2, b1) == -1);
    BlockPtr new_block = Block::join(b1, b2, 1);
    BOOST_CHECK(new_block->size() == 2);
    BOOST_CHECK(new_block->front()->length() == 4);
}

BOOST_AUTO_TEST_CASE (Block_join_bad) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    FragmentPtr f21 = Fragment::create_new(s1, 3, 4);
    FragmentPtr f22 = Fragment::create_new(s2, 3, 4);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    BOOST_CHECK(Block::can_join(b1, b2) == 0);
    BOOST_CHECK(Block::can_join(b2, b1) == 0);
}

BOOST_AUTO_TEST_CASE (Block_try_join) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    FragmentPtr f21 = Fragment::create_new(s1, 3, 4, -1);
    FragmentPtr f22 = Fragment::create_new(s2, 3, 4, -1);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BlockPtr new_block = Block::try_join(b1, b2);
    BOOST_CHECK(new_block);
    BOOST_CHECK(new_block->size() == 2);
    BOOST_CHECK(new_block->front()->length() == 4);
}

BOOST_AUTO_TEST_CASE (Block_try_join_max_gap) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    FragmentPtr f21 = Fragment::create_new(s1, 4, 4, -1);
    FragmentPtr f22 = Fragment::create_new(s2, 4, 4, -1);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BlockPtr new_block = Block::try_join(b1, b2, 0);
    BOOST_CHECK(!new_block);
}

BOOST_AUTO_TEST_CASE (Block_inverse) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b = Block::create_new();
    FragmentPtr f1 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f2 = Fragment::create_new(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    b->inverse();
    BOOST_CHECK(*f1 == Fragment(s1, 1, 2, -1));
    BOOST_CHECK(*f2 == Fragment(s2, 1, 2, -1));
}

BOOST_AUTO_TEST_CASE (Block_patch) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b = Block::create_new();
    Fragment f0(s1, 3, 5, -1);
    FragmentPtr f1 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f2 = Fragment::create_new(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    FragmentDiff diff = f1->diff_to(f0);
    b->patch(diff);
    BOOST_CHECK(*f1 == Fragment(s1, 3, 5, -1));
    BOOST_CHECK(*f2 == Fragment(s2, 3, 5, -1));
}

BOOST_AUTO_TEST_CASE (Block_split) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b = Block::create_new();
    FragmentPtr f1 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f2 = Fragment::create_new(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    BlockPtr new_block = b->split(1);
    BOOST_CHECK(*f1 == Fragment(s1, 1, 1, 1));
    BOOST_CHECK(*f2 == Fragment(s2, 1, 1, 1));
    BOOST_REQUIRE(new_block && new_block->size() == 2);
    BOOST_CHECK(new_block->front()->str() == "g");
    BOOST_CHECK(f1->next());
    BOOST_CHECK(f2->next());
}

BOOST_AUTO_TEST_CASE (Block_max_shift_end) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b = Block::create_new();
    FragmentPtr f1 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f2 = Fragment::create_new(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    BOOST_CHECK(b->max_shift_end() == 15);
    b->inverse();
    BOOST_CHECK(b->max_shift_end() == 1);
}

BOOST_AUTO_TEST_CASE (Block_max_shift_end_two_blocks) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccGAgcggacggcc");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = Fragment::create_new(s1, 11, 12);
    FragmentPtr f22 = Fragment::create_new(s2, 6, 7);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ -1) == 15);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 0) == 3);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 1) == 4);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 5) == 8);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ -1) == 5);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ 0) == 5);
    b1->inverse();
    BOOST_CHECK(b1->max_shift_end(/* overlap */ -1) == 1);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 0) == 1);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ -1) == 5);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ 0) == 5);
    b2->inverse();
    BOOST_CHECK(b1->max_shift_end(/* overlap */ -1) == 1);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 0) == 1);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ -1) == 6);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ 0) == 3);
}

BOOST_AUTO_TEST_CASE (Block_expand_basic) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b = Block::create_new();
    FragmentPtr f1 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f2 = Fragment::create_new(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    b->expand();
    BOOST_CHECK(f1->min_pos() == 0);
    BOOST_CHECK(f1->max_pos() == s1->size() - 1);
    BOOST_CHECK(f2->min_pos() == 0);
    BOOST_CHECK(f2->max_pos() == s2->size() - 1);
}

BOOST_AUTO_TEST_CASE (Block_expand_3) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAaatcagatcg");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAgcggacggcc");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TGGTCCGAgcggacggcc");
    BlockPtr b = Block::create_new();
    FragmentPtr f1 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f2 = Fragment::create_new(s2, 1, 2);
    FragmentPtr f3 = Fragment::create_new(s3, 1, 2);
    b->insert(f1);
    b->insert(f2);
    b->insert(f3);
    PairAligner aliner(1, 0);
    b->expand(&aliner);
    BOOST_CHECK(f1->min_pos() == 0);
    BOOST_CHECK(f1->max_pos() == 7);
    BOOST_CHECK(f2->min_pos() == 0);
    BOOST_CHECK(f2->max_pos() == 7);
    BOOST_CHECK(f3->min_pos() == 0);
    BOOST_CHECK(f3->max_pos() == 7);
}

BOOST_AUTO_TEST_CASE (Block_expand_two_blocks) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccGAcggccgcgga");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = Fragment::create_new(s1, 11, 12);
    FragmentPtr f22 = Fragment::create_new(s2, 6, 7);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    b1->expand(); // overlap = 0
    BOOST_CHECK(f11->min_pos() == 0);
    BOOST_CHECK(f11->max_pos() == 5);
    BOOST_CHECK(f12->min_pos() == 0);
    BOOST_CHECK(f12->max_pos() == 5);
    b2->expand(); // overlap = 0
    BOOST_CHECK(f21->min_pos() == 11);
    BOOST_CHECK(f21->max_pos() == 17);
    BOOST_CHECK(f22->min_pos() == 6);
    BOOST_CHECK(f22->max_pos() == 12);
    PairAligner eq(0);
    b1->expand(&eq, 100, 0, /* max_overlap */ 1);
    BOOST_CHECK(f11->min_pos() == 0);
    BOOST_CHECK(f11->max_pos() == 6);
    BOOST_CHECK(f12->min_pos() == 0);
    BOOST_CHECK(f12->max_pos() == 6);
}

BOOST_AUTO_TEST_CASE (Block_expand_intersection) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGTCCGAGCGGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGTCCGAGCGGAcggcc");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 5);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 5);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = Fragment::create_new(s1, 3, 12);
    FragmentPtr f22 = Fragment::create_new(s2, 3, 12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    b1->expand(); // overlap = 0
    BOOST_CHECK(f11->min_pos() == 0);
    BOOST_CHECK(f11->max_pos() == 5);
    BOOST_CHECK(f12->min_pos() == 0);
    BOOST_CHECK(f12->max_pos() == 5);
    b2->expand(); // overlap = 0
    BOOST_CHECK(f21->min_pos() == 3);
    BOOST_CHECK(f21->max_pos() == 17);
    BOOST_CHECK(f22->min_pos() == 3);
    BOOST_CHECK(f22->max_pos() == 17);
}

BOOST_AUTO_TEST_CASE (Block_common_positions) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgacggccgcgga");
    BlockPtr b1 = Block::create_new();
    b1->insert(Fragment::create_new(s1, 1, 2));
    b1->insert(Fragment::create_new(s2, 1, 2));
    BOOST_CHECK(b1->common_positions(Fragment(s1, 10, 11)) == 0);
    BOOST_CHECK(b1->common_positions(Fragment(s1, 2, 5)) == 1);
}

BOOST_AUTO_TEST_CASE (Block_expand_blocks_by_fragments) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = Fragment::create_new(s1, 11, 12);
    b2->insert(f21);
    Fragment::connect(f11, f21);
    BOOST_CHECK(b2->expand_by_fragments());
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
}

BOOST_AUTO_TEST_CASE (Block_expand_blocks_by_fragments_batch_1) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = Fragment::create_new(s1, 11, 12);
    b2->insert(f21);
    Fragment::connect(f11, f21);
    BOOST_CHECK(b2->expand_by_fragments(/* aligner */ 0, /* batch */ 1));
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
}

BOOST_AUTO_TEST_CASE (Block_expand_blocks_by_fragments_length_1) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGgtccgagcgGacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGgtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 1);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 1);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = Fragment::create_new(s1, 11, 11);
    b2->insert(f21);
    Fragment::connect(f11, f21);
    BOOST_CHECK(b2->expand_by_fragments());
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
}

BOOST_AUTO_TEST_CASE (Block_expand_blocks_by_fragments_high) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcg");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcg");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = Fragment::create_new(s1, 1, 2);
    FragmentPtr f12 = Fragment::create_new(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    const int SEQ_NUMBER = 1000;
    std::vector<SequencePtr> seqs;
    for (int i = 2; i < SEQ_NUMBER; i++) {
        SequencePtr s = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcg");
        seqs.push_back(s);
        b1->insert(Fragment::create_new(s, 1, 2));
    }
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = Fragment::create_new(s1, 11, 12);
    FragmentPtr f22 = Fragment::create_new(s2, 11, 12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BOOST_CHECK(b2->expand_by_fragments());
    BOOST_CHECK(b2->size() == SEQ_NUMBER);
}

BOOST_AUTO_TEST_CASE (Block_expand_blocks_by_fragments_self_neighbour) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GaGaGaGaG");
    FragmentPtr f11 = Fragment::create_new(s1, 0, 0);
    FragmentPtr f12 = Fragment::create_new(s1, 2, 2);
    FragmentPtr f13 = Fragment::create_new(s1, 4, 4);
    FragmentPtr f14 = Fragment::create_new(s1, 6, 6);
    FragmentPtr f15 = Fragment::create_new(s1, 8, 8);
    BlockPtr b = Block::create_new();
    b->insert(f11);
    b->insert(f12);
    b->insert(f13);
    b->insert(f14);
    b->insert(f15);
    std::vector<SequencePtr> seqs;
    for (int i = 0; i < 1000; i++) {
        SequencePtr s = boost::make_shared<InMemorySequence>("gagaGagag");
        seqs.push_back(s);
        b->insert(Fragment::create_new(s, 4, 4));
    }
    Fragment::connect(f11, f12);
    Fragment::connect(f12, f13);
    Fragment::connect(f13, f14);
    Fragment::connect(f14, f15);
    b->expand_by_fragments(); // check that no segfault occurs
}

BOOST_AUTO_TEST_CASE (Block_merge) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GaGaGaGaG");
    FragmentPtr f11 = Fragment::create_new(s1, 0, 0);
    FragmentPtr f12 = Fragment::create_new(s1, 2, 2);
    FragmentPtr f13 = Fragment::create_new(s1, 4, 4);
    FragmentPtr f21 = Fragment::create_new(s1, 4, 4, -1);
    FragmentPtr f22 = Fragment::create_new(s1, 6, 6, -1);
    FragmentPtr f23 = Fragment::create_new(s1, 8, 8, -1);
    BlockPtr b1 = Block::create_new();
    b1->insert(f11);
    b1->insert(f12);
    b1->insert(f13);
    BlockPtr b2 = Block::create_new();
    b2->insert(f21);
    b2->insert(f22);
    b2->insert(f23);
    b1->merge(b2);
    BOOST_CHECK(b2->empty());
    BOOST_CHECK(b1->size() == 5);
    BOOST_CHECK(b1->identity() == 1);
}

