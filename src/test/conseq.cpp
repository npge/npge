/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include "Sequence.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "ConSeq.hpp"
#include "DeConSeq.hpp"

bool near(size_t a, size_t b) {
    return a == b || a - b == 1 || b - a == 1;
}

BOOST_AUTO_TEST_CASE (ConSeq_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<CompactSequence>("-caggccgg");
    SequencePtr s2 = boost::make_shared<CompactSequence>("-caggctg-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("gctggatg-");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 1, s3->size() - 1);
    Block* block = new Block;
    block->insert(f1);
    block->insert(f2);
    block->insert(f3);
    block->set_name("block1");
    BlockSetPtr block_set = new_bs();
    block_set->insert(block);
    ConSeq conseq(block_set);
    conseq.set_empty_block_set();
    conseq.run();
    BOOST_REQUIRE(conseq.block_set()->seqs().size() == 1);
    SequencePtr new_seq = conseq.block_set()->seqs()[0];
    BOOST_CHECK(new_seq->name() == "block1");
    BOOST_CHECK(new_seq->contents() == "caggctgg");
    Fragment* f4 = new Fragment(new_seq, 0, 2); // cag
    Fragment* f5 = new Fragment(new_seq, 4, 6, -1); // cag
    Block* block_l2 = new Block;
    block_l2->insert(f4);
    block_l2->insert(f5);
    BlockSetPtr block_set_l2 = new_bs();
    block_set_l2->insert(block_l2);
    DeConSeq deconseq(block_set_l2);
    deconseq.set_empty_block_set();
    deconseq.run();
    BOOST_REQUIRE(deconseq.block_set()->size() == 1);
    Block* slice = *deconseq.block_set()->begin();
    BOOST_CHECK(slice->size() == 6);
    BOOST_CHECK(slice->alignment_length() == 3);
    BOOST_CHECK(slice->consensus_string() == "cag");
    BOOST_CHECK(slice->identity() < 1);
    BOOST_FOREACH (Fragment* f, *slice) {
        size_t min_pos = f->min_pos();
        size_t max_pos = f->max_pos();
        int ori = f->ori();
        if (f->seq() == s3.get()) {
            min_pos -= 1;
            max_pos -= 1;
        }
        BOOST_CHECK((near(min_pos, 0) && near(max_pos, 2) && ori == 1) ||
                    (near(min_pos, 4) && near(max_pos, 6) && ori == -1));
    }
}

BOOST_AUTO_TEST_CASE (ConSeq_alignment) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<CompactSequence>("-caggccgg");
    SequencePtr s2 = boost::make_shared<CompactSequence>("-caggctg-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("gctggatg-");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 1, s3->size() - 1);
    f1->set_row(new CompactAlignmentRow("caggccgg"));
    f2->set_row(new CompactAlignmentRow("caggctg-"));
    f3->set_row(new CompactAlignmentRow("ctggatg-"));
    Block* block = new Block;
    block->insert(f1);
    block->insert(f2);
    block->insert(f3);
    BlockSetPtr block_set = new_bs();
    block_set->insert(block);
    ConSeq conseq(block_set);
    conseq.set_empty_block_set();
    conseq.run();
    BOOST_REQUIRE(conseq.block_set()->seqs().size() == 1);
    SequencePtr new_seq = conseq.block_set()->seqs()[0];
    BOOST_CHECK(new_seq->contents() == "caggctgg");
    Fragment* f4 = new Fragment(new_seq, 0, 2); // cag
    Fragment* f5 = new Fragment(new_seq, 4, 6, -1); // cag
    Block* block_l2 = new Block;
    block_l2->insert(f4);
    block_l2->insert(f5);
    BlockSetPtr block_set_l2 = new_bs();
    block_set_l2->insert(block_l2);
    DeConSeq deconseq(block_set_l2);
    deconseq.set_empty_block_set();
    deconseq.run();
    BOOST_REQUIRE(deconseq.block_set()->size() == 1);
    Block* slice = *deconseq.block_set()->begin();
    BOOST_FOREACH (Fragment* f, *slice) {
        size_t min_pos = f->min_pos();
        size_t max_pos = f->max_pos();
        int ori = f->ori();
        if (f->seq() == s3.get()) {
            min_pos -= 1;
            max_pos -= 1;
        }
        BOOST_CHECK((min_pos == 0 && max_pos == 2 && ori == 1) ||
                    (min_pos == 4 && max_pos == 6 && ori == -1));
    }
}

