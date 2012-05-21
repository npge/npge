/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"

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

