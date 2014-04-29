/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
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
    SequencePtr s1 = boost::make_shared<CompactSequence>("-CAGGCCGG");
    SequencePtr s2 = boost::make_shared<CompactSequence>("-CAGGCTG-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("GCTGGATG-");
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
    BOOST_CHECK(new_seq->contents() == "CAGGCCGG");
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
    BOOST_CHECK(slice->consensus_string() == "CAG");
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
    SequencePtr s1 = boost::make_shared<CompactSequence>("-CAGGCCGG");
    SequencePtr s2 = boost::make_shared<CompactSequence>("-CAGGCTG-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("GCTGGATG-");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 1, s3->size() - 1);
    f1->set_row(new CompactAlignmentRow("CAGGCCGG"));
    f2->set_row(new CompactAlignmentRow("CAGGCTG-"));
    f3->set_row(new CompactAlignmentRow("CTGGATG-"));
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
    BOOST_CHECK(new_seq->contents() == "CAGGCTGG");
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

BOOST_AUTO_TEST_CASE (DeConSeq_alignment) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<CompactSequence>("-CAGGCCGG");
    SequencePtr s2 = boost::make_shared<CompactSequence>("-CAGGCTG-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("GCTGGATG-");
    s1->set_name("s1");
    s2->set_name("s2");
    s3->set_name("s3");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    f1->set_row(new CompactAlignmentRow("-CAGGCCGG"));
    f2->set_row(new CompactAlignmentRow("-CAGGCTG-"));
    f3->set_row(new CompactAlignmentRow("GCTGGATG-"));
    Block* b = new Block;
    b->insert(f1);
    b->insert(f2);
    b->insert(f3);
    SequencePtr s4((new CompactSequence("GCAGAGCCGG")));
    s4->set_name("s4");
    Fragment* f4 = new Fragment(s4, 0, s4->size() - 1);
    f4->set_row(new CompactAlignmentRow("GCAGAGCCGG"));
    Block* ba = new Block;
    ba->insert(f4);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b);
    block_set->insert(ba);
    ConSeq conseq(block_set);
    conseq.set_empty_block_set();
    conseq.run();
    BOOST_REQUIRE(conseq.block_set()->seqs().size() == 2);
    SequencePtr new_seq = conseq.block_set()->seqs()[0];
    SequencePtr new_seq_a = conseq.block_set()->seqs()[1];
    if (new_seq->contents() == "GCAGAGCCGG") {
        new_seq.swap(new_seq_a);
    }
    BOOST_REQUIRE(new_seq->contents() == "GCAGGCTGG");
    BOOST_REQUIRE(new_seq_a->contents() == "GCAGAGCCGG");
    Fragment* cf1 = new Fragment(new_seq, 0, new_seq->size() - 1);
    Fragment* cf2 = new Fragment(new_seq_a, 0, new_seq_a->size() - 1);
    cf1->set_row(new CompactAlignmentRow("GCAG-GCTGG"));
    cf2->set_row(new CompactAlignmentRow("GCAGAGCCGG"));
    Block* cb = new Block;
    cb->insert(cf1);
    cb->insert(cf2);
    BlockSetPtr cbs = new_bs();
    cbs->insert(cb);
    DeConSeq deconseq(cbs);
    deconseq.set_empty_block_set();
    deconseq.run();
    BlockSetPtr result = deconseq.block_set();
    BOOST_REQUIRE(result->size() == 1);
    Block* result_block = result->front();
    BOOST_REQUIRE(result_block->size() == 4);
    Strings aln;
    BOOST_FOREACH (Fragment* f, *result_block) {
        aln.push_back(f->str());
    }
    std::sort(aln.begin(), aln.end());
    BOOST_CHECK(aln[0] == "-CAG-GCCGG");
    BOOST_CHECK(aln[1] == "-CAG-GCTG-");
    BOOST_CHECK(aln[2] == "GCAGAGCCGG");
    BOOST_CHECK(aln[3] == "GCTG-GATG-");
}

