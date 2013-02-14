/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "BlocksExpander.hpp"
#include "Connector.hpp"

BOOST_AUTO_TEST_CASE (BlocksExpander_expand_block_by_fragments) {
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
    Fragment::connect(f11, f21);
    BlocksExpander expander;
    BOOST_CHECK(expander.expand(b2));
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (BlocksExpander_expand_block_by_fragments_batch_1) {
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
    Fragment::connect(f11, f21);
    BlocksExpander expander(/* batch */ 1);
    BOOST_CHECK(expander.expand(b2));
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (BlocksExpander_expand_block_by_fragments_length_1) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGgtccgagcgGacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGgtccgagcggacggcc");
    Block* b1 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 1);
    Fragment* f12 = new Fragment(s2, 1, 1);
    b1->insert(f11);
    b1->insert(f12);
    Block* b2 = new Block();
    Fragment* f21 = new Fragment(s1, 11, 11);
    b2->insert(f21);
    Fragment::connect(f11, f21);
    BlocksExpander expander;
    BOOST_CHECK(expander.expand(b2));
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (BlocksExpander_expand_block_by_fragments_high) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcg");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcg");
    Block* b1 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 2);
    Fragment* f12 = new Fragment(s2, 1, 2);
    b1->insert(f11);
    b1->insert(f12);
    const int SEQ_NUMBER = 1000;
    std::vector<SequencePtr> seqs;
    for (int i = 2; i < SEQ_NUMBER; i++) {
        SequencePtr s = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcg");
        seqs.push_back(s);
        b1->insert(new Fragment(s, 1, 2));
    }
    Block* b2 = new Block();
    Fragment* f21 = new Fragment(s1, 11, 12);
    Fragment* f22 = new Fragment(s2, 11, 12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    BlocksExpander expander;
    BOOST_CHECK(expander.expand(b2));
    BOOST_CHECK(b2->size() == SEQ_NUMBER);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (BlocksExpander_expand_block_by_fragments_self_neighbor) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GaGaGaGaG");
    Fragment* f11 = new Fragment(s1, 0, 0);
    Fragment* f12 = new Fragment(s1, 2, 2);
    Fragment* f13 = new Fragment(s1, 4, 4);
    Fragment* f14 = new Fragment(s1, 6, 6);
    Fragment* f15 = new Fragment(s1, 8, 8);
    Block* b = new Block();
    b->insert(f11);
    b->insert(f12);
    b->insert(f13);
    b->insert(f14);
    b->insert(f15);
    std::vector<SequencePtr> seqs;
    for (int i = 0; i < 1000; i++) {
        SequencePtr s = boost::make_shared<InMemorySequence>("gagaGagag");
        seqs.push_back(s);
        b->insert(new Fragment(s, 4, 4));
    }
    Fragment::connect(f11, f12);
    Fragment::connect(f12, f13);
    Fragment::connect(f13, f14);
    Fragment::connect(f14, f15);
    BlocksExpander expander;
    expander.expand(b); // check that no segfault occurs
    delete b;
}

BOOST_AUTO_TEST_CASE (BlocksExpander_expand_blocks_by_fragments) {
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
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    BlocksExpander expander;
    BOOST_CHECK(expander.apply(block_set));
    BOOST_CHECK(!expander.expand(b2));
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
}

BOOST_AUTO_TEST_CASE (BlocksExpander_expand_blocks_by_fragments_batch_1) {
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
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    BlocksExpander expander(/* batch */ 1);
    BOOST_CHECK(expander.apply(block_set));
    BOOST_CHECK(!expander.expand(b2));
    BOOST_CHECK(b2->size() == 2);
    BOOST_CHECK(f12->next());
}

