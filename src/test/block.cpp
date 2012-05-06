/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
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
    BOOST_REQUIRE(block->front() == f1 || block->front() == f2);
    BOOST_REQUIRE(f1->block() == block);
    BOOST_REQUIRE(f2->block() == block);
    BOOST_FOREACH (FragmentPtr fragment, *block) {
        BOOST_REQUIRE(fragment == f1 || fragment == f2);
    }
    block->erase(f1);
    BOOST_REQUIRE(block->size() == 1);
    BOOST_REQUIRE(!f1->block());
    BOOST_REQUIRE(f2->block());
    BOOST_REQUIRE(!block->has(f1));
    BOOST_REQUIRE(block->has(f2));
    BOOST_REQUIRE(block->front() == f2);
}

