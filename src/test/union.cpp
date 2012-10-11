/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Union.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

BOOST_AUTO_TEST_CASE (Union_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    Block* b = new Block();
    b->insert(f1);
    b->insert(f2);
    BlockSetPtr source = boost::make_shared<BlockSet>();
    source->insert(b);
    Union u(source);
    BlockSetPtr dest = boost::make_shared<BlockSet>();
    u.apply(dest);
    BOOST_REQUIRE(dest->size() == 1);
    BOOST_CHECK(dest->front()->size() == 2);
}

