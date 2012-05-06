/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"

BOOST_AUTO_TEST_CASE (Block_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccgagcggacggcc");
    FragmentPtr f1 = boost::make_shared<Fragment>(s1, 2, 6, 1);
    FragmentPtr f2 = boost::make_shared<Fragment>(s1, 9, 13, -1);
    BlockPtr block = boost::make_shared<Block>();
    block->insert(f1);
    block->insert(f2);
    BOOST_REQUIRE(block->size() == 2);
    BOOST_REQUIRE(f1->str() == f2->str());
}

