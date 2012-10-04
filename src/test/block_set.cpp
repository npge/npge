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
#include "BlockSet.hpp"
#include "Joiner.hpp"
#include "Filter.hpp"
#include "Connector.hpp"

BOOST_AUTO_TEST_CASE (BlockSet_connect) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Block* b1 = new Block();
    Block* b2 = new Block();
    Block* b3 = new Block();
    b1->insert(f1);
    b2->insert(f2);
    b3->insert(f3);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Connector connector;
    connector.apply(block_set);
    BOOST_CHECK(f1->next() == f2);
    BOOST_CHECK(f2->prev() == f1);
    BOOST_CHECK(f2->next() == f3);
    BOOST_CHECK(f3->prev() == f2);
}

BOOST_AUTO_TEST_CASE (BlockSet_clone) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Block* b1 = new Block();
    Block* b2 = new Block();
    Block* b3 = new Block();
    b1->insert(f1);
    b2->insert(f2);
    b3->insert(f3);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Connector connector;
    connector.apply(block_set);
    BlockSetPtr block_set_copy = block_set->clone();
    BOOST_CHECK(block_set_copy->size() == 3);
    BOOST_CHECK(block_set_copy->front()->front()->prev() ||
                block_set_copy->front()->front()->next());
}

BOOST_AUTO_TEST_CASE (BlockSet_filter) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 4, 6, -1);
    Fragment* f3 = new Fragment(s1, 7, 8, 1);
    Block* b1 = new Block();
    Block* b2 = new Block();
    Block* b3 = new Block();
    b1->insert(f1);
    b2->insert(f2);
    b3->insert(f3);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Filter filter(3, 1);
    filter.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
}

BOOST_AUTO_TEST_CASE (BlockSet_expand) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccGAcggccgcgga");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("tGGtccgacggccgcgga");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    b1->insert(new Fragment(s3, 1, 2));
    Block* b2 = new Block();
    b2->insert(new Fragment(s1, 11, 12));
    b2->insert(new Fragment(s2, 6, 7));
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    block_set->expand_blocks();
    BOOST_CHECK(b1->front()->length() == 6);
    BOOST_CHECK(b1->front()->min_pos() == 0);
    BOOST_CHECK(b1->front()->str() == "tggtcc");
    BOOST_CHECK(b2->front()->length() == 7);
    BOOST_CHECK(b2->front()->str() == "gacggcc");
}

BOOST_AUTO_TEST_CASE (BlockSet_expand_blocks_by_fragments) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b1 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    Block* b2 = new Block();
    Fragment* f21 = new Fragment(s1, 11, 12);
    b2->insert(f21);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    BOOST_CHECK(block_set->expand_blocks_by_fragments());
    BOOST_CHECK(!b2->expand_by_fragments());
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
}

BOOST_AUTO_TEST_CASE (BlockSet_expand_blocks_by_fragments_batch_1) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b1 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    Block* b2 = new Block();
    Fragment* f21 = new Fragment(s1, 11, 12);
    b2->insert(f21);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    BOOST_CHECK(block_set->expand_blocks_by_fragments(0, /* batch */ 1));
    BOOST_CHECK(!b2->expand_by_fragments());
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
}

BOOST_AUTO_TEST_CASE (BlockSet_rest) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b1 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    Block* b2 = new Block();
    Fragment* f21 = new Fragment(s1, 11, 12);
    b2->insert(f21);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    BlockSetPtr rest = block_set->rest();
    BOOST_CHECK(rest->size() == 5);
    Filter filter;
    filter.set_min_block_size(1);
    filter.set_min_fragment_length(2);
    filter.apply(rest);
    BOOST_CHECK(rest->size() == 3);
    filter.set_min_fragment_length(6);
    filter.apply(rest);
    BOOST_CHECK(rest->size() == 2);
    filter.set_min_fragment_length(8);
    filter.apply(rest);
    BOOST_CHECK(rest->size() == 2);
    filter.set_min_fragment_length(9);
    filter.apply(rest);
    BOOST_CHECK(rest->size() == 1);
}

