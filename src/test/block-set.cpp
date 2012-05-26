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
    BOOST_CHECK(f1->next() == f2);
    BOOST_CHECK(f2->prev() == f1);
    BOOST_CHECK(f2->next() == f3);
    BOOST_CHECK(f3->prev() == f2);
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
    BOOST_CHECK(block_set->size() == 1);
}

BOOST_AUTO_TEST_CASE (BlockSet_merge) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 1, 2, 1);
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 4, 6, -1);
    FragmentPtr f31 = boost::make_shared<Fragment>(s1, 7, 8, 1);
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 1, 2, 1);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 4, 6, -1);
    FragmentPtr f32 = boost::make_shared<Fragment>(s2, 7, 8, 1);
    BlockPtr b1 = boost::make_shared<Block>();
    BlockPtr b2 = boost::make_shared<Block>();
    BlockPtr b3 = boost::make_shared<Block>();
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
    block_set->connect_fragments();
    block_set->merge();
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->front()->length() == 8);
}

BOOST_AUTO_TEST_CASE (BlockSet_expand) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccGAcggccgcgga");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("tGGtccgacggccgcgga");
    BlockPtr b1 = Block::create_new();
    b1->insert(boost::make_shared<Fragment>(s1, 1, 2));
    b1->insert(boost::make_shared<Fragment>(s2, 1, 2));
    b1->insert(boost::make_shared<Fragment>(s3, 1, 2));
    BlockPtr b2 = Block::create_new();
    b2->insert(boost::make_shared<Fragment>(s1, 11, 12));
    b2->insert(boost::make_shared<Fragment>(s2, 6, 7));
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->connect_fragments();
    block_set->expand_blocks();
    BOOST_CHECK(b1->front()->length() == 6);
    BOOST_CHECK(b1->front()->min_pos() == 0);
    BOOST_CHECK(b1->front()->str() == "tggtcc");
    BOOST_CHECK(b2->front()->length() == 7);
    BOOST_CHECK(b2->front()->str() == "gacggcc");
}

BOOST_AUTO_TEST_CASE (BlockSet_intersections) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ctgc|ACGC|gacgt");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("ctgc|ACGCGA|cgt");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("ctgcac|GCGA|cgt");
    SequencePtr s4 = boost::make_shared<InMemorySequence>("ctgcac|GCGA|cgt");
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 4, 7, -1);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 4, 7, -1);
    BlockPtr b1 = Block::create_new();
    b1->insert(f11);
    b1->insert(f12);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 6, 9);
    FragmentPtr f23 = boost::make_shared<Fragment>(s3, 6, 9);
    FragmentPtr f24 = boost::make_shared<Fragment>(s4, 6, 9);
    BlockPtr b2 = Block::create_new();
    b2->insert(f22);
    b2->insert(f23);
    b2->insert(f24);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->connect_fragments();
    BOOST_CHECK(block_set->intersections());
    f22->set_min_pos(8);
    BOOST_CHECK(!block_set->intersections());
}

BOOST_AUTO_TEST_CASE (BlockSet_resolve_intersections_2) {
    using namespace bloomrepeats;
    /*
    Input:
        Block 1:
            seq1: ---xxxx----
            seq2: ---xxxx----
        Block 2
            seq2: -----xxxx--
            seq3: -----xxxx--
            seq4: -----xxxx--
    Output of resolve_intersections(2):
        Block 1:
            seq1: ---xx------
            seq2: ---xx------
        Block 2
            seq2: -------xx--
            seq3: -------xx--
            seq4: -------xx--
        Block 3:
            seq1: -----xx----
            seq2: -----xx----
            seq3: -----xx----
            seq4: -----xx----
    */
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ctgc|ACAG|gacgt");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    SequencePtr s4 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 4, 7, -1);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 4, 7, -1);
    BlockPtr b1 = Block::create_new();
    b1->insert(f11);
    b1->insert(f12);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 6, 9);
    FragmentPtr f23 = boost::make_shared<Fragment>(s3, 6, 9);
    FragmentPtr f24 = boost::make_shared<Fragment>(s4, 6, 9);
    BlockPtr b2 = Block::create_new();
    b2->insert(f22);
    b2->insert(f23);
    b2->insert(f24);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->connect_fragments();
    block_set->resolve_intersections(2);
    BOOST_REQUIRE(block_set->size() == 3);
    bool b[5] = {0, 0, 0, 0, 0};
    BOOST_FOREACH (BlockPtr block, *block_set) {
        BOOST_REQUIRE(block->size() <= 4);
        b[block->size()] = true;
        if (block->size() == 2) {
            BOOST_CHECK(block->front()->str() == "ac" ||
                        block->front()->str() == "gt");
        } else if (block->size() == 3) {
            BOOST_CHECK(block->front()->str() == "ga" ||
                        block->front()->str() == "tc");
        } else if (block->size() == 4) {
            BOOST_CHECK(block->front()->str() == "ag" ||
                        block->front()->str() == "ct");
        } else {
            BOOST_ERROR("Bad block size");
        }
    }
    BOOST_CHECK(!b[0] && !b[1] && b[2] && b[3] && b[4]);
    BOOST_CHECK(f12->next()->next() == f22);
}

BOOST_AUTO_TEST_CASE (BlockSet_resolve_intersections_3) {
    using namespace bloomrepeats;
    /*
    Input:
        Block 1:
            seq1: ---xxxx----
            seq2: ---xxxx----
        Block 2
            seq2: -----xxxx--
            seq3: -----xxxx--
            seq4: -----xxxx--
    Output of resolve_intersections(3):
        Block 1:
            seq1: ---xx------
            seq2: ---xx------
        Block 2
            seq2: -----xxxx--
            seq3: -----xxxx--
            seq4: -----xxxx--
    */
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ctgc|ACAG|gacgt");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    SequencePtr s4 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 4, 7, -1);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 4, 7, -1);
    BlockPtr b1 = Block::create_new();
    b1->insert(f11);
    b1->insert(f12);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 6, 9);
    FragmentPtr f23 = boost::make_shared<Fragment>(s3, 6, 9);
    FragmentPtr f24 = boost::make_shared<Fragment>(s4, 6, 9);
    BlockPtr b2 = Block::create_new();
    b2->insert(f22);
    b2->insert(f23);
    b2->insert(f24);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->connect_fragments();
    block_set->resolve_intersections(3);
    BOOST_REQUIRE(block_set->size() == 2);
    bool b[4] = {0, 0, 0, 0};
    BOOST_FOREACH (BlockPtr block, *block_set) {
        BOOST_REQUIRE(block->size() <= 3);
        b[block->size()] = true;
        if (block->size() == 2) {
            BOOST_CHECK(block->front()->str() == "ac" ||
                        block->front()->str() == "gt");
        } else if (block->size() == 3) {
            BOOST_CHECK(block->front()->str() == "agga" ||
                        block->front()->str() == "tcct");
        } else {
            BOOST_ERROR("Bad block size");
        }
    }
    BOOST_CHECK(!b[0] && !b[1] && b[2] && b[3]);
    BOOST_CHECK(f12->next() == f22);
}

BOOST_AUTO_TEST_CASE (BlockSet_resolve_intersections_internal) {
    using namespace bloomrepeats;
    /*
    Input:
        Block 1:
            seq1: ---xxxx----
            seq2: ---xxxx----
        Block 2
            seq1: -----xxxx--
            seq2: -----xxxx--
            seq3: -----xxxx--
            seq4: -----xxxx--
    Output of resolve_intersections(2):
        Block 1:
            seq1: ---xx------
            seq2: ---xx------
        Block 2
            seq1: -----xxxx--
            seq2: -----xxxx--
            seq3: -----xxxx--
            seq4: -----xxxx--
    */
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    SequencePtr s4 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    FragmentPtr f11 = boost::make_shared<Fragment>(s1, 4, 7, -1);
    FragmentPtr f12 = boost::make_shared<Fragment>(s2, 4, 7, -1);
    BlockPtr b1 = Block::create_new();
    b1->insert(f11);
    b1->insert(f12);
    FragmentPtr f21 = boost::make_shared<Fragment>(s1, 6, 9);
    FragmentPtr f22 = boost::make_shared<Fragment>(s2, 6, 9);
    FragmentPtr f23 = boost::make_shared<Fragment>(s3, 6, 9);
    FragmentPtr f24 = boost::make_shared<Fragment>(s4, 6, 9);
    BlockPtr b2 = Block::create_new();
    b2->insert(f21);
    b2->insert(f22);
    b2->insert(f23);
    b2->insert(f24);
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->connect_fragments();
    block_set->resolve_intersections(2);
    BOOST_REQUIRE(block_set->size() == 2);
    bool b[5] = {0, 0, 0, 0, 0};
    BOOST_FOREACH (BlockPtr block, *block_set) {
        BOOST_REQUIRE(block->size() <= 4);
        b[block->size()] = true;
        if (block->size() == 2) {
            BOOST_CHECK(block->front()->str() == "ac" ||
                        block->front()->str() == "gt");
        } else if (block->size() == 4) {
            BOOST_CHECK(block->front()->str() == "agga" ||
                        block->front()->str() == "tcct");
        } else {
            BOOST_ERROR("Bad block size");
        }
    }
    BOOST_CHECK(!b[0] && !b[1] && b[2] && !b[3] && b[4]);
    BOOST_CHECK(f11->next() == f21);
    BOOST_CHECK(f12->next() == f22);
}

