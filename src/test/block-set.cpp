/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

BOOST_AUTO_TEST_CASE (BlockSet_connect) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 1, 2, 1);
    FragmentPtr f2 = boost::make_shared<Fragment>(s1, 5, 6, -1);
    FragmentPtr f3 = boost::make_shared<Fragment>(s1, 7, 8, 1);
    BlockPtr b1 = boost::make_shared<Block>();
    BlockPtr b2 = boost::make_shared<Block>();
    BlockPtr b3 = boost::make_shared<Block>();
    b1->insert(f1);
    b2->insert(f2);
    b3->insert(f3);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    block_set->connect_fragments();
    BOOST_REQUIRE(f1->next() == f2);
    BOOST_REQUIRE(f2->prev() == f1);
    BOOST_REQUIRE(f2->next() == f3);
    BOOST_REQUIRE(f3->prev() == f2);
}

BOOST_AUTO_TEST_CASE (BlockSet_filter) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 1, 2, 1);
    FragmentPtr f2 = boost::make_shared<Fragment>(s1, 4, 6, -1);
    FragmentPtr f3 = boost::make_shared<Fragment>(s1, 7, 8, 1);
    BlockPtr b1 = boost::make_shared<Block>();
    BlockPtr b2 = boost::make_shared<Block>();
    BlockPtr b3 = boost::make_shared<Block>();
    b1->insert(f1);
    b2->insert(f2);
    b3->insert(f3);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    block_set->filter(3, 1);
    BOOST_REQUIRE(block_set->size() == 1);
}

