/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Joiner.hpp"
#include "Connector.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

BOOST_AUTO_TEST_CASE (Joiner_Fragment_join) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, 1);
    Fragment* f3 = new Fragment(s1, 7, 8, -1);
    Fragment::connect(f1, f2);
    Fragment::connect(f2, f3);
    Fragment::connect(f3, f1);
    BOOST_CHECK(Joiner::can_join(f1, f2));
    BOOST_CHECK(Joiner::can_join(f2, f1));
    BOOST_CHECK(!Joiner::can_join(f1, f3));
    BOOST_CHECK(!Joiner::can_join(f2, f3));
    Fragment* f12 = Joiner::join(f1, f2);
    BOOST_CHECK(f12->ori() == 1);
    BOOST_CHECK(f12->seq() == s1.get());
    BOOST_CHECK(f12->min_pos() == 1);
    BOOST_CHECK(f12->max_pos() == 6);
    BOOST_CHECK(f12->is_neighbor(*f3));
    BOOST_CHECK(!f12->is_neighbor(*f1));
    BOOST_CHECK(!f12->is_neighbor(*f2));
    delete f1;
    delete f2;
    delete f3;
    delete f12;
}

BOOST_AUTO_TEST_CASE (Joiner_Block_join) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    Block* b1 = new Block();
    Block* b2 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    Fragment* f21 = new Fragment(s1, 3, 4);
    Fragment* f22 = new Fragment(s2, 3, 4);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BOOST_CHECK(Joiner::can_join(b1, b2) == 1);
    BOOST_CHECK(Joiner::can_join(b2, b1) == -1);
    Block* new_block = Joiner::join(b1, b2, 1);
    BOOST_CHECK(new_block->size() == 2);
    BOOST_CHECK(new_block->front()->length() == 4);
    delete b1;
    delete b2;
    delete new_block;
}

BOOST_AUTO_TEST_CASE (Joiner_Block_join_bad) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    Block* b1 = new Block();
    Block* b2 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    Fragment* f21 = new Fragment(s1, 3, 4);
    Fragment* f22 = new Fragment(s2, 3, 4);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    BOOST_CHECK(Joiner::can_join(b1, b2) == 0);
    BOOST_CHECK(Joiner::can_join(b2, b1) == 0);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Joiner_Block_try_join) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    Block* b1 = new Block();
    Block* b2 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    Fragment* f21 = new Fragment(s1, 3, 4, -1);
    Fragment* f22 = new Fragment(s2, 3, 4, -1);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    Joiner joiner;
    Block* new_block = joiner.try_join(b1, b2);
    BOOST_CHECK(new_block);
    BOOST_CHECK(new_block->size() == 2);
    BOOST_CHECK(new_block->front()->length() == 4);
    delete b1;
    delete b2;
    delete new_block;
}

BOOST_AUTO_TEST_CASE (Joiner_fragment) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 1, 2);
    Fragment f2(s1, 5, 6);
    Fragment f3(s1, 8, 8);
    Fragment::connect(&f1, &f2);
    Fragment::connect(&f2, &f3);
    Joiner always_true;
    BOOST_CHECK(always_true.can_join_fragments(&f1, &f2));
    BOOST_CHECK(always_true.can_join_fragments(&f2, &f3));
    Joiner dist_1(1);
    BOOST_CHECK(!dist_1.can_join_fragments(&f1, &f2));
    BOOST_CHECK(dist_1.can_join_fragments(&f2, &f3));
    Joiner ratio_1(-1, 1);
    BOOST_CHECK(ratio_1.can_join_fragments(&f1, &f2));
    BOOST_CHECK(ratio_1.can_join_fragments(&f2, &f3));
    Joiner ratio_05(-1, 0.5);
    BOOST_CHECK(!ratio_05.can_join_fragments(&f1, &f2));
    BOOST_CHECK(!ratio_05.can_join_fragments(&f2, &f3));
}

BOOST_AUTO_TEST_CASE (Joiner_block) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgagatgcgggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tg-gatgcgggcc");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 0, 0));
    b1->insert(new Fragment(s2, 0, 0));
    Block* b2 = new Block();
    b2->insert(new Fragment(s1, 3, 4));
    b2->insert(new Fragment(s2, 2, 3));
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    Joiner always_true;
    BOOST_CHECK(always_true.can_join_blocks(b1, b2));
    Joiner dist_1(1);
    BOOST_CHECK(!dist_1.can_join_blocks(b1, b2));
    Joiner dist_2(2);
    BOOST_CHECK(dist_2.can_join_blocks(b1, b2));
    Joiner gap_ratio_1(-1, -1, 1);
    BOOST_CHECK(!gap_ratio_1.can_join_blocks(b1, b2));
    Joiner gap_ratio_2(-1, -1, 2);
    BOOST_CHECK(gap_ratio_2.can_join_blocks(b1, b2));
    BOOST_CHECK(!gap_ratio_1.can_join_blocks(b1, b2));
}

BOOST_AUTO_TEST_CASE (Joiner_Block_try_join_max_gap) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    Block* b1 = new Block();
    Block* b2 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    Fragment* f21 = new Fragment(s1, 4, 4, -1);
    Fragment* f22 = new Fragment(s2, 4, 4, -1);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    Joiner dist_0(0);
    Block* new_block = dist_0.try_join(b1, b2);
    BOOST_CHECK(!new_block);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Joiner_BlockSet_join) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f11 = new Fragment(s1, 1, 2, 1);
    Fragment* f21 = new Fragment(s1, 4, 6, -1);
    Fragment* f31 = new Fragment(s1, 7, 8, 1);
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f12 = new Fragment(s2, 1, 2, 1);
    Fragment* f22 = new Fragment(s2, 4, 6, -1);
    Fragment* f32 = new Fragment(s2, 7, 8, 1);
    Block* b1 = new Block();
    Block* b2 = new Block();
    Block* b3 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    b3->insert(f31);
    b3->insert(f32);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Connector connector;
    connector.apply(block_set);
    Joiner joiner;
    joiner.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->front()->length() == 8);
}

BOOST_AUTO_TEST_CASE (Joiner_BlockSet_join_max_gap) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f11 = new Fragment(s1, 1, 2, 1);
    Fragment* f21 = new Fragment(s1, 4, 6, -1);
    Fragment* f31 = new Fragment(s1, 7, 8, 1);
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f12 = new Fragment(s2, 1, 2, 1);
    Fragment* f22 = new Fragment(s2, 4, 6, -1);
    Fragment* f32 = new Fragment(s2, 7, 8, 1);
    Block* b1 = new Block();
    Block* b2 = new Block();
    Block* b3 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    b3->insert(f31);
    b3->insert(f32);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Connector connector;
    connector.apply(block_set);
    Joiner dist_0(0);
    dist_0.apply(block_set);
    BOOST_CHECK(block_set->size() == 2);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->front()->length() == 5 ||
                block_set->front()->front()->length() == 2);
}

BOOST_AUTO_TEST_CASE (Joiner_BlockSet_join_wrong) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f11 = new Fragment(s1, 1, 2, 1);
    Fragment* f21 = new Fragment(s1, 4, 6, -1);
    Fragment* f31 = new Fragment(s1, 7, 8, 1);
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f12 = new Fragment(s2, 1, 2, 1);
    Fragment* f22 = new Fragment(s2, 4, 6, -1);
    Fragment* f32 = new Fragment(s2, 7, 8, -1);
    Block* b1 = new Block();
    Block* b2 = new Block();
    Block* b3 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    b3->insert(f31);
    b3->insert(f32);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Connector connector;
    connector.apply(block_set);
    Joiner joiner;
    joiner.apply(block_set);
    BOOST_CHECK(block_set->size() == 2);
}

