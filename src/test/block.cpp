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
#include "Joiner.hpp"

BOOST_AUTO_TEST_CASE (Block_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    Fragment* f1 = new Fragment(s1, 2, 6, 1);
    Fragment* f2 = new Fragment(s1, 9, 13, -1);
    BOOST_REQUIRE(!f1->block());
    Block* block = new Block();
    BOOST_REQUIRE(!block->front());
    block->insert(f1);
    block->insert(f2);
    BOOST_REQUIRE(block->size() == 2);
    BOOST_REQUIRE(block->has(f1));
    BOOST_REQUIRE(block->has(f2));
    BOOST_CHECK(block->front() == f1 || block->front() == f2);
    BOOST_REQUIRE(f1->block() == block);
    BOOST_REQUIRE(f2->block() == block);
    BOOST_FOREACH (Fragment* fragment, *block) {
        BOOST_CHECK(fragment == f1 || fragment == f2);
    }
    block->erase(f1);
    BOOST_CHECK(block->size() == 1);
    BOOST_CHECK(f2->block());
    BOOST_CHECK(!block->has(f1));
    BOOST_CHECK(block->has(f2));
    BOOST_CHECK(block->front() == f2);
    {
        Block* block_copy = block->clone();
        BOOST_CHECK(block_copy->size() == 1);
        BOOST_CHECK(f2->block());
        BOOST_CHECK(!block_copy->has(f2));
        BOOST_CHECK(block_copy->front() != f2);
        BOOST_CHECK(*block_copy->front() == *f2);
        delete block_copy;
    }
    block->clear();
    BOOST_CHECK(block->size() == 0);
    BOOST_CHECK(block->empty());
    delete block;
}

BOOST_AUTO_TEST_CASE (Block_identity) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtTcgagcAgTcggcc");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 0, s1->size() - 1));
    b1->insert(new Fragment(s2, 0, s1->size() - 1));
    BOOST_CHECK(std::abs(b1->identity() - 15. / 18.) < 0.01);
    Block* b2 = new Block();
    b2->insert(new Fragment(s1, 0, 2));
    b2->insert(new Fragment(s2, 0, 5));
    BOOST_CHECK(std::abs(b2->identity() - 0.5) < 0.01);
    Block* b3 = new Block();
    b3->insert(new Fragment(s1, 0, 2));
    b3->insert(new Fragment(s2, 0, 2));
    BOOST_CHECK(std::abs(b3->identity() - 1) < 0.01);
    delete b1;
    delete b2;
    delete b3;
}

BOOST_AUTO_TEST_CASE (Block_match) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    Block* b1 = new Block();
    Block* b2 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    b2->insert(new Fragment(s1, 5, 6));
    b2->insert(new Fragment(s2, 5, 6));
    BOOST_CHECK(Block::match(b1, b2) == 1);
    b1->insert(new Fragment(s1, 8, 9));
    BOOST_CHECK(!Block::match(b1, b2));
    b2->insert(new Fragment(s1, 10, 11));
    BOOST_CHECK(Block::match(b1, b2) == 1);
    b1->insert(new Fragment(s1, 12, 13));
    b2->insert(new Fragment(s1, 14, 15, -1));
    BOOST_CHECK(!Block::match(b1, b2));
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Block_match_1) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    Block* b1 = new Block();
    Block* b2 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    b2->insert(new Fragment(s1, 5, 6, -1));
    b2->insert(new Fragment(s2, 5, 6, -1));
    BOOST_CHECK(Block::match(b1, b2) == -1);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Block_inverse) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    b->inverse();
    BOOST_CHECK(*f1 == Fragment(s1, 1, 2, -1));
    BOOST_CHECK(*f2 == Fragment(s2, 1, 2, -1));
    delete b;
}

BOOST_AUTO_TEST_CASE (Block_patch) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment f0(s1, 3, 5, -1);
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    FragmentDiff diff = f1->diff_to(f0);
    b->patch(diff);
    BOOST_CHECK(*f1 == Fragment(s1, 3, 5, -1));
    BOOST_CHECK(*f2 == Fragment(s2, 3, 5, -1));
    delete b;
}

BOOST_AUTO_TEST_CASE (Block_split) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    Block* new_block = b->split(1);
    BOOST_CHECK(*f1 == Fragment(s1, 1, 1, 1));
    BOOST_CHECK(*f2 == Fragment(s2, 1, 1, 1));
    BOOST_REQUIRE(new_block && new_block->size() == 2);
    BOOST_CHECK(new_block->front()->str() == "g");
    BOOST_CHECK(f1->next());
    BOOST_CHECK(f2->next());
    delete b;
    delete new_block;
}

BOOST_AUTO_TEST_CASE (Block_max_shift_end) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    BOOST_CHECK(b->max_shift_end() == 15);
    b->inverse();
    BOOST_CHECK(b->max_shift_end() == 1);
    delete b;
}

BOOST_AUTO_TEST_CASE (Block_max_shift_end_two_blocks) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccGAgcggacggcc");
    Block* b1 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    Block* b2 = new Block();
    Fragment* f21 = new Fragment(s1, 11, 12);
    Fragment* f22 = new Fragment(s2, 6, 7);
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
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Block_common_positions) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgacggccgcgga");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    BOOST_CHECK(b1->common_positions(Fragment(s1, 10, 11)) == 0);
    BOOST_CHECK(b1->common_positions(Fragment(s1, 2, 5)) == 1);
    delete b1;
}

BOOST_AUTO_TEST_CASE (Block_merge) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GaGaGaGaG");
    Fragment* f11 = new Fragment(s1, 0, 0);
    Fragment* f12 = new Fragment(s1, 2, 2);
    Fragment* f13 = new Fragment(s1, 4, 4);
    Fragment* f21 = new Fragment(s1, 4, 4, -1);
    Fragment* f22 = new Fragment(s1, 6, 6, -1);
    Fragment* f23 = new Fragment(s1, 8, 8, -1);
    Block* b1 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    b1->insert(f13);
    Block* b2 = new Block();
    b2->insert(f21);
    b2->insert(f22);
    b2->insert(f23);
    b1->merge(b2);
    BOOST_CHECK(b2->empty());
    BOOST_CHECK(b1->size() == 5);
    BOOST_CHECK(b1->identity() == 1);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Block_name) {
    using namespace bloomrepeats;
    Block* b = new Block;
    BOOST_CHECK(!b->name().empty());
    b->set_name("abc");
    BOOST_CHECK(b->name() == "abc");
    b->set_name("eee");
    delete b;
    Block block("name");
    BOOST_CHECK(block.name() == "name");
}

