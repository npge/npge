/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cmath>
#include <vector>
#include <sstream>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Block.hpp"
#include "PairAligner.hpp"
#include "Joiner.hpp"
#include "block_stat.hpp"
#include "char_to_size.hpp"
#include "block_hash.hpp"
#include "global.hpp"

BOOST_AUTO_TEST_CASE (Block_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Fragment* f1 = new Fragment(s1, 2, 6, 1);
    Fragment* f2 = new Fragment(s1, 9, 13, -1);
    BOOST_REQUIRE(!f1->block());
    Block* block = new Block();
    BOOST_REQUIRE(!block->front());
    block->insert(f1);
    block->insert(f2);
    BOOST_REQUIRE(block->size() == 2);
    BOOST_REQUIRE(block->has(f1));
    BOOST_REQUIRE(block->has(f2));
    BOOST_CHECK(block->front() == f1 || block->front() == f2);
    BOOST_REQUIRE(f1->block() == block);
    BOOST_REQUIRE(f2->block() == block);
    BOOST_FOREACH (Fragment* fragment, *block) {
        BOOST_CHECK(fragment == f1 || fragment == f2);
    }
    block->erase(f1);
    BOOST_CHECK(block->size() == 1);
    BOOST_CHECK(f2->block());
    BOOST_CHECK(!block->has(f1));
    BOOST_CHECK(block->has(f2));
    BOOST_CHECK(block->front() == f2);
    block->clear();
    BOOST_CHECK(block->size() == 0);
    BOOST_CHECK(block->empty());
    delete block;
}

BOOST_AUTO_TEST_CASE (Block_length) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGT");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGT-CC--");
    Block* b = new Block();
    b->insert(new Fragment(s1, 0, s1->size() - 1));
    BOOST_CHECK(b->alignment_length() == s1->size());
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    f2->set_row(new MapAlignmentRow("TGGT-CC--"));
    b->insert(f2);
    BOOST_CHECK(f2->row()->length() == 9);
    BOOST_CHECK(b->alignment_length() == f2->row()->length());
    delete b;
}

BOOST_AUTO_TEST_CASE (Block_length2) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<CompactSequence>("CAGGACGG");
    SequencePtr s2 = boost::make_shared<CompactSequence>("CAGGAAG-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("CTGGACG-");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    Block block;
    block.insert(f1);
    block.insert(f2);
    block.insert(f3);
    BOOST_CHECK(block.alignment_length() == std::string("CAGGACGG").size());
    BOOST_CHECK(block.consensus_string() == "CAGGACGG");
}

BOOST_AUTO_TEST_CASE (Block_alignment_stat) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TAGTCCG-");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGTT-CG-");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TG---CG-");
    Block b;
    b.insert(new Fragment(s1, 0, s1->size() - 1));
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TGTT-CG-", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    new MapAlignmentRow("TG---CG-", f3);
    b.insert(f3);
    AlignmentStat stat;
    make_stat(stat, &b);
    BOOST_CHECK(stat.ident_nogap() == 3);
    BOOST_CHECK(stat.ident_gap() == 2);
    BOOST_CHECK(stat.noident_nogap() == 1);
    BOOST_CHECK(stat.noident_gap() == 1);
    BOOST_CHECK(stat.pure_gap() == 1);
    BOOST_CHECK(stat.total() == 8);
    BOOST_CHECK(stat.letter_count('A') == 1);
    BOOST_CHECK(stat.letter_count('T') == 6);
    BOOST_CHECK(stat.letter_count('G') == 6);
    BOOST_CHECK(stat.letter_count('C') == 4);
    BOOST_CHECK(stat.gc() > 0.5);
    //
    AlignmentStat stat5;
    make_stat(stat5, &b, /* start */ 5);
    BOOST_CHECK(stat5.ident_nogap() == 2);
    BOOST_CHECK(stat5.ident_gap() == 0);
    BOOST_CHECK(stat5.noident_nogap() == 0);
    BOOST_CHECK(stat5.noident_gap() == 0);
    BOOST_CHECK(stat5.pure_gap() == 1);
    BOOST_CHECK(stat5.total() == 3);
    BOOST_CHECK(stat5.letter_count('A') == 0);
    BOOST_CHECK(stat5.letter_count('T') == 0);
    BOOST_CHECK(stat5.letter_count('G') == 3);
    BOOST_CHECK(stat5.letter_count('C') == 3);
    BOOST_CHECK(stat5.gc() > 0.99);
    //
    AlignmentStat stat56;
    make_stat(stat56, &b, /* start */ 5, /* stop */ 6);
    BOOST_CHECK(stat56.ident_nogap() == 2);
    BOOST_CHECK(stat56.ident_gap() == 0);
    BOOST_CHECK(stat56.noident_nogap() == 0);
    BOOST_CHECK(stat56.noident_gap() == 0);
    BOOST_CHECK(stat56.pure_gap() == 0);
    BOOST_CHECK(stat56.total() == 2);
    BOOST_CHECK(stat56.letter_count('A') == 0);
    BOOST_CHECK(stat56.letter_count('T') == 0);
    BOOST_CHECK(stat56.letter_count('G') == 3);
    BOOST_CHECK(stat56.letter_count('C') == 3);
    BOOST_CHECK(stat56.gc() > 0.99);
    //
    bool ident, gap, pure_gap;
    int atgc[LETTERS_NUMBER];
    for (int i = 0; i < LETTERS_NUMBER; i++) {
        atgc[i] = 0;
    }
    test_column(&b, 0, ident, gap, pure_gap, atgc);
    BOOST_CHECK(ident);
    BOOST_CHECK(!gap);
    BOOST_CHECK(!pure_gap);
    BOOST_CHECK(atgc[A] == 0);
    BOOST_CHECK(atgc[T] == 3);
    BOOST_CHECK(atgc[G] == 0);
    BOOST_CHECK(atgc[C] == 0);
    BOOST_CHECK(atgc[N] == 0);
    //
    for (int i = 0; i < LETTERS_NUMBER; i++) {
        atgc[i] = 0;
    }
    test_column(&b, 1, ident, gap, pure_gap, atgc);
    BOOST_CHECK(!ident);
    BOOST_CHECK(!gap);
    BOOST_CHECK(!pure_gap);
    BOOST_CHECK(atgc[A] == 1);
    BOOST_CHECK(atgc[T] == 0);
    BOOST_CHECK(atgc[G] == 2);
    BOOST_CHECK(atgc[C] == 0);
    BOOST_CHECK(atgc[N] == 0);
    //
    for (int i = 0; i < LETTERS_NUMBER; i++) {
        atgc[i] = 0;
    }
    test_column(&b, 2, ident, gap, pure_gap, atgc);
    BOOST_CHECK(!ident);
    BOOST_CHECK(gap);
    BOOST_CHECK(!pure_gap);
    BOOST_CHECK(atgc[A] == 0);
    BOOST_CHECK(atgc[T] == 1);
    BOOST_CHECK(atgc[G] == 1);
    BOOST_CHECK(atgc[C] == 0);
    BOOST_CHECK(atgc[N] == 0);
    //
    for (int i = 0; i < LETTERS_NUMBER; i++) {
        atgc[i] = 0;
    }
    test_column(&b, 3, ident, gap, pure_gap, atgc);
    BOOST_CHECK(ident);
    BOOST_CHECK(gap);
    BOOST_CHECK(!pure_gap);
    BOOST_CHECK(atgc[A] == 0);
    BOOST_CHECK(atgc[T] == 2);
    BOOST_CHECK(atgc[G] == 0);
    BOOST_CHECK(atgc[C] == 0);
    BOOST_CHECK(atgc[N] == 0);
    //
    for (int i = 0; i < LETTERS_NUMBER; i++) {
        atgc[i] = 0;
    }
    test_column(&b, 4, ident, gap, pure_gap, atgc);
    BOOST_CHECK(ident);
    BOOST_CHECK(gap);
    BOOST_CHECK(!pure_gap);
    BOOST_CHECK(atgc[A] == 0);
    BOOST_CHECK(atgc[T] == 0);
    BOOST_CHECK(atgc[G] == 0);
    BOOST_CHECK(atgc[C] == 1);
    BOOST_CHECK(atgc[N] == 0);
    //
    for (int i = 0; i < LETTERS_NUMBER; i++) {
        atgc[i] = 0;
    }
    test_column(&b, 7, ident, gap, pure_gap, atgc);
    BOOST_CHECK(ident);
    BOOST_CHECK(gap);
    BOOST_CHECK(pure_gap);
    BOOST_CHECK(atgc[A] == 0);
    BOOST_CHECK(atgc[T] == 0);
    BOOST_CHECK(atgc[G] == 0);
    BOOST_CHECK(atgc[C] == 0);
    BOOST_CHECK(atgc[N] == 0);
}

BOOST_AUTO_TEST_CASE (Block_slice) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TAGTCCG-");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGTT-CG-");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TG---CG-");
    Block b;
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    new MapAlignmentRow("TAGTCCG-", f1);
    b.insert(f1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TGTT-CG-", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    new MapAlignmentRow("TG---CG-", f3);
    b.insert(f3);
    boost::scoped_ptr<Block> b46((b.slice(4, 6)));
    BOOST_REQUIRE(b46->size() == 3);
    Strings ff;
    BOOST_FOREACH (Fragment* f, *b46) {
        ff.push_back(f->str());
    }
    std::sort(ff.begin(), ff.end());
    BOOST_CHECK(ff[0] == "-CG");
    BOOST_CHECK(ff[1] == "-CG");
    BOOST_CHECK(ff[2] == "CCG");
    AlignmentStat stat46;
    make_stat(stat46, b46.get());
    BOOST_CHECK(stat46.ident_nogap() == 2);
    BOOST_CHECK(stat46.ident_gap() == 1);
    BOOST_CHECK(stat46.noident_nogap() == 0);
    BOOST_CHECK(stat46.noident_gap() == 0);
    BOOST_CHECK(stat46.pure_gap() == 0);
    BOOST_CHECK(stat46.total() == 3);
    BOOST_CHECK(stat46.letter_count('A') == 0);
    BOOST_CHECK(stat46.letter_count('T') == 0);
    BOOST_CHECK(stat46.letter_count('G') == 3);
    BOOST_CHECK(stat46.letter_count('C') == 4);
    BOOST_CHECK(stat46.gc() > 0.99);
}

BOOST_AUTO_TEST_CASE (Block_slice_gaps) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TAGTCCG-");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGTT-CG-");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TG---CG-");
    Block b;
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    new MapAlignmentRow("TAGTCCG-", f1);
    b.insert(f1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TGTT-CG-", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    new MapAlignmentRow("TG---CG-", f3);
    b.insert(f3);
    boost::scoped_ptr<Block> b26((b.slice(2, 6)));
    BOOST_REQUIRE(b26->size() == 3);
    Strings ff;
    BOOST_FOREACH (Fragment* f, *b26) {
        ff.push_back(f->str());
    }
    std::sort(ff.begin(), ff.end());
    BOOST_CHECK(ff[0] == "---CG");
    BOOST_CHECK(ff[1] == "GTCCG");
    BOOST_CHECK(ff[2] == "TT-CG");
}

BOOST_AUTO_TEST_CASE (Block_slice_reverse) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TAGTCCG-");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGTT-CG-");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TG---CG-");
    Block b;
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    new MapAlignmentRow("TAGTCCG-", f1);
    b.insert(f1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TGTT-CG-", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    new MapAlignmentRow("TG---CG-", f3);
    b.insert(f3);
    boost::scoped_ptr<Block> b63((b.slice(6, 3)));
    BOOST_REQUIRE(b63->size() == 3);
    Strings ff;
    BOOST_FOREACH (Fragment* f, *b63) {
        ff.push_back(f->str());
    }
    std::sort(ff.begin(), ff.end());
    BOOST_CHECK(ff[0] == "CG--");
    BOOST_CHECK(ff[1] == "CG-A");
    BOOST_CHECK(ff[2] == "CGGA");
    AlignmentStat stat63;
    make_stat(stat63, b63.get());
    BOOST_CHECK(stat63.ident_nogap() == 2);
    BOOST_CHECK(stat63.ident_gap() == 2);
    BOOST_CHECK(stat63.noident_nogap() == 0);
    BOOST_CHECK(stat63.noident_gap() == 0);
    BOOST_CHECK(stat63.pure_gap() == 0);
    BOOST_CHECK(stat63.total() == 4);
    BOOST_CHECK(stat63.letter_count('A') == 2);
    BOOST_CHECK(stat63.letter_count('T') == 0);
    BOOST_CHECK(stat63.letter_count('G') == 4);
    BOOST_CHECK(stat63.letter_count('C') == 3);
    BOOST_CHECK(stat63.gc() > 0.5);
}

BOOST_AUTO_TEST_CASE (Block_weak) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TAGTCCG-");
    Block strong_block;
    strong_block.insert(new Fragment(s1, 0, s1->size() - 1));
    Block weak_block;
    weak_block.set_weak(true);
    weak_block.insert(strong_block.front());
    BOOST_CHECK(!weak_block.empty());
    BOOST_CHECK(weak_block.size() == 1);
    BOOST_CHECK(!strong_block.empty());
    BOOST_CHECK(strong_block.size() == 1);
    BOOST_CHECK(weak_block.front()->block() == &strong_block);
    weak_block.clear();
    BOOST_CHECK(!strong_block.empty());
    BOOST_CHECK(strong_block.size() == 1);
}

BOOST_AUTO_TEST_CASE (Block_weak2) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TAGTCCG-");
    Block strong_block;
    strong_block.insert(new Fragment(s1, 0, s1->size() - 1));
    Block weak_block;
    weak_block.set_weak(true);
    weak_block.insert(strong_block.front());
    BOOST_CHECK(!weak_block.empty());
    BOOST_CHECK(weak_block.size() == 1);
    BOOST_CHECK(!strong_block.empty());
    BOOST_CHECK(strong_block.size() == 1);
    BOOST_CHECK(weak_block.front()->block() == &strong_block);
    weak_block.set_weak(false);
    BOOST_CHECK(!weak_block.empty());
    BOOST_CHECK(weak_block.size() == 1);
    BOOST_CHECK(weak_block.weak() == false);
    BOOST_CHECK(!strong_block.empty());
    BOOST_CHECK(strong_block.size() == 1);
    BOOST_CHECK(strong_block.weak() == true);
    BOOST_CHECK(weak_block.front()->block() == &weak_block);
}

BOOST_AUTO_TEST_CASE (Block_identity) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtTcgagcAgTcggcc");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 0, s1->size() - 1));
    b1->insert(new Fragment(s2, 0, s1->size() - 1));
    BOOST_CHECK(std::abs(b1->identity() - 15. / 18.) < 0.01);
    Block* b2 = new Block();
    b2->insert(new Fragment(s1, 0, 2));
    b2->insert(new Fragment(s2, 0, 5));
    BOOST_CHECK(std::abs(b2->identity() - 1.0) < 0.01);
    Block* b3 = new Block();
    b3->insert(new Fragment(s1, 0, 2));
    b3->insert(new Fragment(s2, 0, 2));
    BOOST_CHECK(std::abs(b3->identity() - 1) < 0.01);
    delete b1;
    delete b2;
    delete b3;
}

BOOST_AUTO_TEST_CASE (Block_consensus) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TAGTCCG-");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGTT-CG-");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TG---CG-");
    Block b;
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    new MapAlignmentRow("TAGTCCG-", f1);
    b.insert(f1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TGTT-CG-", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s2->size() - 1);
    new MapAlignmentRow("TG---CG-", f3);
    b.insert(f3);
    std::stringstream consensus_stream;
    b.consensus(consensus_stream, /* gap */ 'A');
    std::string consensus_string = consensus_stream.str();
    BOOST_CHECK(consensus_string == "TGGTCCGA" ||
                consensus_string == "TGTTCCGA");
    BOOST_CHECK(b.consensus_string() == "TGGTCCGA" ||
                b.consensus_string() == "TGTTCCGA");
}

BOOST_AUTO_TEST_CASE (Block_match) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block* b1 = new Block();
    Block* b2 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    b2->insert(new Fragment(s1, 5, 6));
    b2->insert(new Fragment(s2, 5, 6));
    BOOST_CHECK(Block::match(b1, b2) == 1);
    b1->insert(new Fragment(s1, 8, 9));
    BOOST_CHECK(!Block::match(b1, b2));
    b2->insert(new Fragment(s1, 10, 11));
    BOOST_CHECK(Block::match(b1, b2) == 1);
    b1->insert(new Fragment(s1, 12, 13));
    b2->insert(new Fragment(s1, 14, 15, -1));
    BOOST_CHECK(!Block::match(b1, b2));
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Block_match_1) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block* b1 = new Block();
    Block* b2 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    b2->insert(new Fragment(s1, 5, 6, -1));
    b2->insert(new Fragment(s2, 5, 6, -1));
    BOOST_CHECK(Block::match(b1, b2) == -1);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Block_inverse) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    b->inverse();
    BOOST_CHECK(*f1 == Fragment(s1, 1, 2, -1));
    BOOST_CHECK(*f2 == Fragment(s2, 1, 2, -1));
    delete b;
}

BOOST_AUTO_TEST_CASE (Block_patch) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment f0(s1, 3, 5, -1);
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    FragmentDiff diff = f1->diff_to(f0);
    b->patch(diff);
    BOOST_CHECK(*f1 == Fragment(s1, 3, 5, -1));
    BOOST_CHECK(*f2 == Fragment(s2, 3, 5, -1));
    delete b;
}

BOOST_AUTO_TEST_CASE (Block_split) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    Block* new_block = b->split(1);
    BOOST_CHECK(*f1 == Fragment(s1, 1, 1, 1));
    BOOST_CHECK(*f2 == Fragment(s2, 1, 1, 1));
    BOOST_REQUIRE(new_block && new_block->size() == 2);
    BOOST_CHECK(new_block->front()->str() == "G");
    BOOST_CHECK(f1->next());
    BOOST_CHECK(f2->next());
    delete b;
    delete new_block;
}

BOOST_AUTO_TEST_CASE (Block_max_shift_end) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    BOOST_CHECK(b->max_shift_end() == 15);
    b->inverse();
    BOOST_CHECK(b->max_shift_end() == 1);
    delete b;
}

BOOST_AUTO_TEST_CASE (Block_max_shift_end_two_blocks) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccGAgcggacggcc");
    Block* b1 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    Block* b2 = new Block();
    Fragment* f21 = new Fragment(s1, 11, 12);
    Fragment* f22 = new Fragment(s2, 6, 7);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ -1) == 15);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 0) == 3);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 1) == 4);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 5) == 8);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ -1) == 5);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ 0) == 5);
    b1->inverse();
    BOOST_CHECK(b1->max_shift_end(/* overlap */ -1) == 1);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 0) == 1);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ -1) == 5);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ 0) == 5);
    b2->inverse();
    BOOST_CHECK(b1->max_shift_end(/* overlap */ -1) == 1);
    BOOST_CHECK(b1->max_shift_end(/* overlap */ 0) == 1);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ -1) == 6);
    BOOST_CHECK(b2->max_shift_end(/* overlap */ 0) == 3);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Block_common_positions) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgacggccgcgga");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    BOOST_CHECK(b1->common_positions(Fragment(s1, 10, 11)) == 0);
    BOOST_CHECK(b1->common_positions(Fragment(s1, 2, 5)) == 1);
    delete b1;
}

BOOST_AUTO_TEST_CASE (Block_merge) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GaGaGaGaG");
    Fragment* f11 = new Fragment(s1, 0, 0);
    Fragment* f12 = new Fragment(s1, 2, 2);
    Fragment* f13 = new Fragment(s1, 4, 4);
    Fragment* f21 = new Fragment(s1, 4, 4, -1);
    Fragment* f22 = new Fragment(s1, 6, 6, -1);
    Fragment* f23 = new Fragment(s1, 8, 8, -1);
    Block* b1 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    b1->insert(f13);
    Block* b2 = new Block();
    b2->insert(f21);
    b2->insert(f22);
    b2->insert(f23);
    b1->merge(b2);
    BOOST_CHECK(b2->empty());
    BOOST_CHECK(b1->size() == 5);
    BOOST_CHECK(b1->identity() == 1);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Block_name) {
    using namespace bloomrepeats;
    Block* b = new Block;
    BOOST_CHECK(!b->name().empty());
    b->set_name("abc");
    BOOST_CHECK(b->name() == "abc");
    b->set_name("eee");
    delete b;
    Block block("name");
    BOOST_CHECK(block.name() == "name");
}

BOOST_AUTO_TEST_CASE (Block_hash) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GaGaGaGaG");
    Fragment* f11 = new Fragment(s1, 0, 4);
    Fragment* f12 = new Fragment(s1, 4, 5, -1);
    boost::scoped_ptr<Block> b1((new Block));
    b1->insert(f11);
    b1->insert(f12);
    boost::scoped_ptr<Block> b2((b1->clone()));
    b2->inverse();
    BOOST_CHECK(block_hash(b1.get()) == block_hash(b2.get()));
}

