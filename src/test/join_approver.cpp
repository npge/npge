/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "JoinApprover.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

BOOST_AUTO_TEST_CASE (JoinApprover_fragment) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagatgcgggcc");
    Fragment f1(s1, 1, 2);
    Fragment f2(s1, 5, 6);
    Fragment f3(s1, 8, 8);
    Fragment::connect(&f1, &f2);
    Fragment::connect(&f2, &f3);
    JoinApprover always_true;
    BOOST_CHECK(always_true.can_join_fragments(&f1, &f2));
    BOOST_CHECK(always_true.can_join_fragments(&f2, &f3));
    JoinApprover dist_1(1);
    BOOST_CHECK(!dist_1.can_join_fragments(&f1, &f2));
    BOOST_CHECK(dist_1.can_join_fragments(&f2, &f3));
}

BOOST_AUTO_TEST_CASE (JoinApprover_block) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgagatgcgggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tg-gatgcgggcc");
    BlockPtr b1 = Block::create_new();
    b1->insert(Fragment::create_new(s1, 0, 0));
    b1->insert(Fragment::create_new(s2, 0, 0));
    BlockPtr b2 = Block::create_new();
    b2->insert(Fragment::create_new(s1, 3, 4));
    b2->insert(Fragment::create_new(s2, 2, 3));
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->connect_fragments();
    JoinApprover always_true;
    BOOST_CHECK(always_true.can_join_blocks(b1, b2));
    JoinApprover dist_1(1);
    BOOST_CHECK(!dist_1.can_join_blocks(b1, b2));
    JoinApprover dist_2(2);
    BOOST_CHECK(dist_2.can_join_blocks(b1, b2));
}

