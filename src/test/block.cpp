/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "PairAligner.hpp"

BOOST_AUTO_TEST_CASE (Block_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 2, 6, 1);
    FragmentPtr f2 = boost::make_shared<Fragment>(s1, 9, 13, -1);
    BOOST_REQUIRE(!f1->block());
    BlockPtr block = boost::make_shared<Block>();
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
    BOOST_CHECK(!f1->block());
    BOOST_CHECK(f2->block());
    BOOST_CHECK(!block->has(f1));
    BOOST_CHECK(block->has(f2));
    BOOST_CHECK(block->front() == f2);
    block->clear();
    BOOST_CHECK(block->size() == 0);
    BOOST_CHECK(block->empty());
    BOOST_CHECK(!f2->block());
}

BOOST_AUTO_TEST_CASE (Block_identity) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtTcgagcAgTcggcc");
    BlockPtr b1 = Block::create_new();
    b1->insert(boost::make_shared<Fragment>(s1, 0, s1->size() - 1));
    b1->insert(boost::make_shared<Fragment>(s2, 0, s1->size() - 1));
    BOOST_CHECK(std::abs(b1->identity() - 15. / 18.) < 0.01);
    BlockPtr b2 = Block::create_new();
    b2->insert(boost::make_shared<Fragment>(s1, 0, 2));
    b2->insert(boost::make_shared<Fragment>(s2, 0, 5));
    BOOST_CHECK(std::abs(b2->identity() - 0.5) < 0.01);
    BlockPtr b3 = Block::create_new();
    b3->insert(boost::make_shared<Fragment>(s1, 0, 2));
    b3->insert(boost::make_shared<Fragment>(s2, 0, 2));
    BOOST_CHECK(std::abs(b3->identity() - 1) < 0.01);
}

BOOST_AUTO_TEST_CASE (Block_match) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    b1->insert(boost::make_shared<Fragment>(s1, 1, 2));
    b1->insert(boost::make_shared<Fragment>(s2, 1, 2));
    b2->insert(boost::make_shared<Fragment>(s1, 5, 6));
    b2->insert(boost::make_shared<Fragment>(s2, 5, 6));
    BOOST_CHECK(Block::match(b1, b2) == 1);
    b1->insert(boost::make_shared<Fragment>(s1, 8, 9));
    BOOST_CHECK(!Block::match(b1, b2));
    b2->insert(boost::make_shared<Fragment>(s1, 10, 11));
    BOOST_CHECK(Block::match(b1, b2) == 1);
    b1->insert(boost::make_shared<Fragment>(s1, 12, 13));
    b2->insert(boost::make_shared<Fragment>(s1, 14, 15, -1));
    BOOST_CHECK(!Block::match(b1, b2));
}

BOOST_AUTO_TEST_CASE (Block_match_1) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    b1->insert(boost::make_shared<Fragment>(s1, 1, 2));
    b1->insert(boost::make_shared<Fragment>(s2, 1, 2));
    b2->insert(boost::make_shared<Fragment>(s1, 5, 6, -1));
    b2->insert(boost::make_shared<Fragment>(s2, 5, 6, -1));
    BOOST_CHECK(Block::match(b1, b2) == -1);
}

BOOST_AUTO_TEST_CASE (Block_merge) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 1, 2);
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 3, 4);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 3, 4);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BOOST_CHECK(Block::can_merge(b1, b2) == 1);
    BOOST_CHECK(Block::can_merge(b2, b1) == -1);
    BlockPtr new_block = Block::merge(b1, b2, 1);
    BOOST_CHECK(new_block->size() == 2);
    BOOST_CHECK(new_block->front()->length() == 4);
}

BOOST_AUTO_TEST_CASE (Block_merge_bad) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 1, 2);
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 3, 4);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 3, 4);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    BOOST_CHECK(Block::can_merge(b1, b2) == 0);
    BOOST_CHECK(Block::can_merge(b2, b1) == 0);
}

BOOST_AUTO_TEST_CASE (Block_try_merge) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    BlockPtr b2 = Block::create_new();
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 1, 2);
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 3, 4, -1);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 3, 4, -1);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BlockPtr new_block = Block::try_merge(b1, b2);
    BOOST_CHECK(new_block);
    BOOST_CHECK(new_block->size() == 2);
    BOOST_CHECK(new_block->front()->length() == 4);
}

BOOST_AUTO_TEST_CASE (Block_inverse) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b = Block::create_new();
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f2 = boost::make_shared<Fragment>(s2, 1, 2);
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
    FragmentPtr f0 = boost::make_shared<Fragment>(s1, 3, 5, -1);
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f2 = boost::make_shared<Fragment>(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    FragmentDiff diff = f1->diff_to(*f0);
    b->patch(diff);
    BOOST_CHECK(*f1 == Fragment(s1, 3, 5, -1));
    BOOST_CHECK(*f2 == Fragment(s2, 3, 5, -1));
}

BOOST_AUTO_TEST_CASE (Block_max_shift_end) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b = Block::create_new();
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f2 = boost::make_shared<Fragment>(s2, 1, 2);
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
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 11, 12);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 6, 7);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ true) == 15);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ false) == 3);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ true) == 5);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ false) == 5);
    b1->inverse();
    BOOST_CHECK(b1->max_shift_end(/* overlap */ true) == 1);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ false) == 1);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ true) == 5);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ false) == 5);
    b2->inverse();
    BOOST_CHECK(b1->max_shift_end(/* overlap */ true) == 1);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ false) == 1);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ true) == 6);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ false) == 3);
}

BOOST_AUTO_TEST_CASE (Block_expand_basic) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b = Block::create_new();
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f2 = boost::make_shared<Fragment>(s2, 1, 2);
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
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f2 = boost::make_shared<Fragment>(s2, 1, 2);
    FragmentPtr f3 = boost::make_shared<Fragment>(s3, 1, 2);
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
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 11, 12);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 6, 7);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    b1->expand(); // overlap = false
    BOOST_CHECK(f11->min_pos() == 0);
    BOOST_CHECK(f11->max_pos() == 5);
    BOOST_CHECK(f12->min_pos() == 0);
    BOOST_CHECK(f12->max_pos() == 5);
    b2->expand(); // overlap = false
    BOOST_CHECK(f21->min_pos() == 11);
    BOOST_CHECK(f21->max_pos() == 17);
    BOOST_CHECK(f22->min_pos() == 6);
    BOOST_CHECK(f22->max_pos() == 12);
}

BOOST_AUTO_TEST_CASE (Block_expand_intersection) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGTCCGAGCGGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGTCCGAGCGGAcggcc");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 1, 5);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 1, 5);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 3, 12);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 3, 12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    b1->expand(); // overlap = false
    BOOST_CHECK(f11->min_pos() == 0);
    BOOST_CHECK(f11->max_pos() == 5);
    BOOST_CHECK(f12->min_pos() == 0);
    BOOST_CHECK(f12->max_pos() == 5);
    b2->expand(); // overlap = false
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
    b1->insert(boost::make_shared<Fragment>(s1, 1, 2));
    b1->insert(boost::make_shared<Fragment>(s2, 1, 2));
    BOOST_CHECK(b1->common_positions(Fragment(s1, 10, 11)) == 0);
    BOOST_CHECK(b1->common_positions(Fragment(s1, 2, 5)) == 1);
}

BOOST_AUTO_TEST_CASE (Block_expand_blocks_by_fragments) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    BlockPtr b1 = Block::create_new();
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 1, 2);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    BlockPtr b2 = Block::create_new();
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 11, 12);
    b2->insert(f21);
    Fragment::connect(f11, f21);
    BOOST_CHECK(b2->expand_by_fragments());
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
}

