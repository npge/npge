/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Joiner.hpp"
#include "Connector.hpp"
#include "OriByMajority.hpp"
#include "Sequence.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

BOOST_AUTO_TEST_CASE (Joiner_Fragment_join) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f1 = new Fragment(s1, 1, 2, 1);
    Fragment* f2 = new Fragment(s1, 5, 6, 1);
    Fragment* f3 = new Fragment(s1, 7, 8, -1);
    Fragment::connect(f1, f2);
    Fragment::connect(f2, f3);
    Fragment::connect(f3, f1);
    BOOST_CHECK(Joiner::can_join(f1, f2));
    BOOST_CHECK(Joiner::can_join(f2, f1));
    BOOST_CHECK(!Joiner::can_join(f1, f3));
    BOOST_CHECK(!Joiner::can_join(f2, f3));
    Fragment* f12 = Joiner::join(f1, f2);
    BOOST_CHECK(f12->ori() == 1);
    BOOST_CHECK(f12->seq() == s1.get());
    BOOST_CHECK(f12->min_pos() == 1);
    BOOST_CHECK(f12->max_pos() == 6);
    BOOST_CHECK(f12->is_neighbor(*f3));
    BOOST_CHECK(!f12->is_neighbor(*f1));
    BOOST_CHECK(!f12->is_neighbor(*f2));
    delete f1;
    delete f2;
    delete f3;
    delete f12;
}

BOOST_AUTO_TEST_CASE (Joiner_Block_join) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block* b1 = new Block();
    Block* b2 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    Fragment* f21 = new Fragment(s1, 3, 4);
    Fragment* f22 = new Fragment(s2, 3, 4);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BOOST_CHECK(Joiner::can_join(b1, b2) == 1);
    BOOST_CHECK(Joiner::can_join(b2, b1) == -1);
    Block* new_block = Joiner().join_blocks(b1, b2, 1);
    BOOST_CHECK(new_block->size() == 2);
    BOOST_CHECK(new_block->front()->length() == 4);
    delete b1;
    delete b2;
    delete new_block;
}

BOOST_AUTO_TEST_CASE (Joiner_Block_join_bad) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block* b1 = new Block();
    Block* b2 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    Fragment* f21 = new Fragment(s1, 3, 4);
    Fragment* f22 = new Fragment(s2, 3, 4);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    BOOST_CHECK(Joiner::can_join(b1, b2) == 0);
    BOOST_CHECK(Joiner::can_join(b2, b1) == 0);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Joiner_Block_try_join) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block* b1 = new Block();
    Block* b2 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    Fragment* f21 = new Fragment(s1, 3, 4, -1);
    Fragment* f22 = new Fragment(s2, 3, 4, -1);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    Joiner joiner;
    Block* new_block = joiner.try_join(b1, b2);
    BOOST_CHECK(new_block);
    BOOST_CHECK(new_block->size() == 2);
    BOOST_CHECK(new_block->front()->length() == 4);
    delete b1;
    delete b2;
    delete new_block;
}

BOOST_AUTO_TEST_CASE (Joiner_fragment) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGATGCGGGCC");
    Fragment f1(s1, 1, 2);
    Fragment f2(s1, 5, 6);
    Fragment f3(s1, 8, 8);
    Fragment::connect(&f1, &f2);
    Fragment::connect(&f2, &f3);
    Joiner always_true;
    BOOST_CHECK(always_true.can_join_fragments(&f1, &f2));
    BOOST_CHECK(always_true.can_join_fragments(&f2, &f3));
}

BOOST_AUTO_TEST_CASE (Joiner_block) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGAGATGCGGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TG-GATGCGGGCC");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 0, 0));
    b1->insert(new Fragment(s2, 0, 0));
    Block* b2 = new Block();
    b2->insert(new Fragment(s1, 3, 4));
    b2->insert(new Fragment(s2, 2, 3));
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    Joiner always_true;
    BOOST_CHECK(always_true.can_join_blocks(b1, b2));
}

BOOST_AUTO_TEST_CASE (Joiner_Block_try_join_max_gap) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block* b1 = new Block();
    Block* b2 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    Fragment* f21 = new Fragment(s1, 4, 4, -1);
    Fragment* f22 = new Fragment(s2, 4, 4, -1);
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (Joiner_BlockSet_join) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f11 = new Fragment(s1, 1, 2, 1);
    Fragment* f21 = new Fragment(s1, 4, 6, -1);
    Fragment* f31 = new Fragment(s1, 7, 8, 1);
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f12 = new Fragment(s2, 1, 2, 1);
    Fragment* f22 = new Fragment(s2, 4, 6, -1);
    Fragment* f32 = new Fragment(s2, 7, 8, 1);
    Block* b1 = new Block();
    Block* b2 = new Block();
    Block* b3 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    b3->insert(f31);
    b3->insert(f32);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Connector connector;
    connector.apply(block_set);
    Joiner joiner;
    joiner.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->front()->length() == 8);
}

BOOST_AUTO_TEST_CASE (Joiner_BlockSet_join_wrong) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f11 = new Fragment(s1, 1, 2, 1);
    Fragment* f21 = new Fragment(s1, 4, 6, -1);
    Fragment* f31 = new Fragment(s1, 7, 8, 1);
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tggtcCGAGATgcgggcc");
    Fragment* f12 = new Fragment(s2, 1, 2, 1);
    Fragment* f22 = new Fragment(s2, 4, 6, -1);
    Fragment* f32 = new Fragment(s2, 7, 8, -1);
    Block* b1 = new Block();
    Block* b2 = new Block();
    Block* b3 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    b3->insert(f31);
    b3->insert(f32);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Connector connector;
    connector.apply(block_set);
    Joiner joiner;
    joiner.apply(block_set);
    BOOST_CHECK(block_set->size() == 2);
}

BOOST_AUTO_TEST_CASE (Joiner_alternation) {
    using namespace npge;
    std::string str = "ATATATATATAT";
    SequencePtr s = boost::make_shared<InMemorySequence>(str);
    Block* b1 = new Block;
    b1->insert(new Fragment(s, 0, 1));
    b1->insert(new Fragment(s, 4, 5));
    Block* b2 = new Block;
    b2->insert(new Fragment(s, 2, 3));
    b2->insert(new Fragment(s, 6, 7));
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    Joiner always_true;
    BOOST_CHECK(always_true.can_join_blocks(b1, b2));
    always_true.apply(block_set);
    BOOST_REQUIRE(block_set->size() == 1);
    Block* b = block_set->front();
    BOOST_CHECK(b->size() == 2);
    BOOST_CHECK(b->alignment_length() == 4);
    BOOST_CHECK(b->consensus_string() == "ATAT");
}

BOOST_AUTO_TEST_CASE (Joiner_alignment) {
    using namespace npge;
    SequencePtr s1((new InMemorySequence("ACTGG")));
    Fragment* f11 = new Fragment(s1, 0, 2, 1);
    f11->set_row(new CompactAlignmentRow("ACT"));
    Fragment* f21 = new Fragment(s1, 3, 4, 1);
    f21->set_row(new CompactAlignmentRow("GG"));
    SequencePtr s2((new InMemorySequence("A-TGG")));
    Fragment* f12 = new Fragment(s2, 0, 1, 1);
    f12->set_row(new CompactAlignmentRow("A-T"));
    Fragment* f22 = new Fragment(s2, 2, 3, 1);
    f22->set_row(new CompactAlignmentRow("GG"));
    Block* b1 = new Block;
    Block* b2 = new Block;
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Joiner joiner;
    joiner.apply(block_set);
    OriByMajority obm;
    obm.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->consensus_string() == "ACTGG");
}

BOOST_AUTO_TEST_CASE (Joiner_alignment_2) {
    using namespace npge;
    SequencePtr s1((new InMemorySequence("ACTGG")));
    Fragment* f11 = new Fragment(s1, 0, 1, 1);
    f11->set_row(new CompactAlignmentRow("AC"));
    Fragment* f21 = new Fragment(s1, 2, 4, 1);
    f21->set_row(new CompactAlignmentRow("TGG"));
    SequencePtr s2((new InMemorySequence("A-TGG")));
    Fragment* f12 = new Fragment(s2, 0, 0, 1);
    f12->set_row(new CompactAlignmentRow("A-"));
    Fragment* f22 = new Fragment(s2, 1, 3, 1);
    f22->set_row(new CompactAlignmentRow("TGG"));
    Block* b1 = new Block;
    Block* b2 = new Block;
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Joiner joiner;
    joiner.apply(block_set);
    OriByMajority obm;
    obm.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->consensus_string() == "ACTGG");
}

BOOST_AUTO_TEST_CASE (Joiner_alignment_rev) {
    using namespace npge;
    SequencePtr s1((new InMemorySequence("ACTGG")));
    Fragment* f11 = new Fragment(s1, 0, 1, -1);
    f11->set_row(new CompactAlignmentRow("GT"));
    Fragment* f21 = new Fragment(s1, 2, 4, 1);
    f21->set_row(new CompactAlignmentRow("TGG"));
    SequencePtr s2((new InMemorySequence("A-TGG")));
    Fragment* f12 = new Fragment(s2, 0, 0, -1);
    f12->set_row(new CompactAlignmentRow("-T"));
    Fragment* f22 = new Fragment(s2, 1, 3, 1);
    f22->set_row(new CompactAlignmentRow("TGG"));
    Block* b1 = new Block;
    Block* b2 = new Block;
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Joiner joiner;
    joiner.apply(block_set);
    OriByMajority obm;
    obm.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->consensus_string() == "ACTGG");
}

BOOST_AUTO_TEST_CASE (Joiner_alignment_rev_2) {
    using namespace npge;
    SequencePtr s1((new InMemorySequence("ACTGG")));
    Fragment* f11 = new Fragment(s1, 0, 1, 1);
    f11->set_row(new CompactAlignmentRow("AC"));
    Fragment* f21 = new Fragment(s1, 2, 4, -1);
    f21->set_row(new CompactAlignmentRow("CCA"));
    SequencePtr s2((new InMemorySequence("A-TGG")));
    Fragment* f12 = new Fragment(s2, 0, 0, 1);
    f12->set_row(new CompactAlignmentRow("A-"));
    Fragment* f22 = new Fragment(s2, 1, 3, -1);
    f22->set_row(new CompactAlignmentRow("CCA"));
    Block* b1 = new Block;
    Block* b2 = new Block;
    b1->insert(f11);
    b1->insert(f12);
    b2->insert(f21);
    b2->insert(f22);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Joiner joiner;
    joiner.apply(block_set);
    OriByMajority obm;
    obm.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->consensus_string() == "ACTGG");
}

using namespace npge;

struct AlignerData {
    SequencePtr s1, s2;
    Fragment* f11, *f21, *f12, *f22;
    Block* b1, *b2;
    BlockSetPtr block_set;

    AlignerData() {
        // ##$$$##
        // 0123456
        // ACTGAAT
        // ACT-AAT
        // 012-345
        // ##$$$##
        s1.reset(new InMemorySequence("ACTGAAT"));
        f11 = new Fragment(s1, 0, 1, 1);
        f11->set_row(new CompactAlignmentRow("AC"));
        f21 = new Fragment(s1, 5, 6, 1);
        f21->set_row(new CompactAlignmentRow("AT"));
        s2.reset(new InMemorySequence("ACT-AAT"));
        f12 = new Fragment(s2, 0, 1, 1);
        f12->set_row(new CompactAlignmentRow("AC"));
        f22 = new Fragment(s2, 4, 5, 1);
        f22->set_row(new CompactAlignmentRow("AT"));
        b1 = new Block;
        b2 = new Block;
        b1->insert(f11);
        b1->insert(f12);
        b2->insert(f21);
        b2->insert(f22);
        block_set = new_bs();
        block_set->insert(b1);
        block_set->insert(b2);
    }
};

BOOST_AUTO_TEST_CASE (Joiner_aligner) {
    using namespace npge;
    AlignerData ad;
    BlockSetPtr block_set = ad.block_set;
    Joiner joiner;
    joiner.apply(block_set);
    OriByMajority obm;
    obm.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->consensus_string() == "ACTGAAT");
}

BOOST_AUTO_TEST_CASE (Joiner_aligner_2) {
    using namespace npge;
    AlignerData ad;
    ad.b1->inverse();
    BlockSetPtr block_set = ad.block_set;
    Joiner joiner;
    joiner.apply(block_set);
    OriByMajority obm;
    obm.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->consensus_string() == "ACTGAAT");
}

BOOST_AUTO_TEST_CASE (Joiner_aligner_3) {
    using namespace npge;
    AlignerData ad;
    ad.b2->inverse();
    BlockSetPtr block_set = ad.block_set;
    Joiner joiner;
    joiner.apply(block_set);
    OriByMajority obm;
    obm.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->consensus_string() == "ACTGAAT");
}

BOOST_AUTO_TEST_CASE (Joiner_aligner_4) {
    using namespace npge;
    AlignerData ad;
    ad.b1->inverse();
    ad.b2->inverse();
    BlockSetPtr block_set = ad.block_set;
    Joiner joiner;
    joiner.apply(block_set);
    OriByMajority obm;
    obm.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->consensus_string() == "ACTGAAT");
}

