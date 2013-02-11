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
#include "ConSeq.hpp"

BOOST_AUTO_TEST_CASE (ConSeq_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<CompactSequence>("caggacgg");
    SequencePtr s2 = boost::make_shared<CompactSequence>("caggaag-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("ctggacg-");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    Block* block = new Block;
    block->insert(f1);
    block->insert(f2);
    block->insert(f3);
    block->set_name("block1");
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->insert(block);
    ConSeq conseq(block_set);
    conseq.set_empty_block_set();
    conseq.run();
    BOOST_REQUIRE(conseq.block_set()->seqs().size() == 1);
    SequencePtr new_seq = conseq.block_set()->seqs()[0];
    BOOST_CHECK(new_seq->name() == "block1");
    BOOST_CHECK(new_seq->contents() == "caggacgg");
}

