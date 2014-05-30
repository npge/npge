/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "OverlapsResolver.hpp"
#include "Connector.hpp"

BOOST_AUTO_TEST_CASE (OverlapsResolver_overlaps) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ctgc|ACGC|gacgt");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("ctgc|ACGCGA|cgt");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("ctgcac|GCGA|cgt");
    SequencePtr s4 = boost::make_shared<InMemorySequence>("ctgcac|GCGA|cgt");
    Fragment* f11 = new Fragment(s1, 4, 7, -1);
    Fragment* f12 = new Fragment(s2, 4, 7, -1);
    Block* b1 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    Fragment* f22 = new Fragment(s2, 6, 9);
    Fragment* f23 = new Fragment(s3, 6, 9);
    Fragment* f24 = new Fragment(s4, 6, 9);
    Block* b2 = new Block();
    b2->insert(f22);
    b2->insert(f23);
    b2->insert(f24);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    OverlapsResolver resolver;
    resolver.set_block_set(block_set);
    BOOST_CHECK(resolver.overlaps());
    f22->set_min_pos(8);
    BOOST_CHECK(!resolver.overlaps());
}

BOOST_AUTO_TEST_CASE (OverlapsResolver_main) {
    using namespace npge;
    /*
    Input:
        Block 1:
            seq1: ---xxxx----
            seq2: ---xxxx----
        Block 2
            seq2: -----xxxx--
            seq3: -----xxxx--
            seq4: -----xxxx--
    Output of OverlapsResolver:
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
    Fragment* f11 = new Fragment(s1, 4, 7, -1);
    Fragment* f12 = new Fragment(s2, 4, 7, -1);
    Block* b1 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    Fragment* f22 = new Fragment(s2, 6, 9);
    Fragment* f23 = new Fragment(s3, 6, 9);
    Fragment* f24 = new Fragment(s4, 6, 9);
    Block* b2 = new Block();
    b2->insert(f22);
    b2->insert(f23);
    b2->insert(f24);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    OverlapsResolver resolver;
    resolver.apply(block_set);
    BOOST_REQUIRE(block_set->size() == 3);
    bool b[5] = {0, 0, 0, 0, 0};
    BOOST_FOREACH (Block* block, *block_set) {
        BOOST_REQUIRE(block->size() <= 4);
        b[block->size()] = true;
        if (block->size() == 2) {
            BOOST_CHECK(block->front()->str() == "AC" ||
                        block->front()->str() == "GT");
        } else if (block->size() == 3) {
            BOOST_CHECK(block->front()->str() == "GA" ||
                        block->front()->str() == "TC");
        } else if (block->size() == 4) {
            BOOST_CHECK(block->front()->str() == "AG" ||
                        block->front()->str() == "CT");
        } else {
            BOOST_ERROR("Bad block size");
        }
    }
    BOOST_CHECK(!b[0] && !b[1] && b[2] && b[3] && b[4]);
}

BOOST_AUTO_TEST_CASE (OverlapsResolver_internal) {
    using namespace npge;
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
    */
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    SequencePtr s4 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    Fragment* f11 = new Fragment(s1, 4, 7, -1);
    Fragment* f12 = new Fragment(s2, 4, 7, -1);
    Block* b1 = new Block();
    b1->insert(f11);
    b1->insert(f12);
    Fragment* f21 = new Fragment(s1, 6, 9);
    Fragment* f22 = new Fragment(s2, 6, 9);
    Fragment* f23 = new Fragment(s3, 6, 9);
    Fragment* f24 = new Fragment(s4, 6, 9);
    Block* b2 = new Block();
    b2->insert(f21);
    b2->insert(f22);
    b2->insert(f23);
    b2->insert(f24);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    OverlapsResolver resolver;
    resolver.apply(block_set);
    BOOST_REQUIRE(block_set->size() >= 2);
}

BOOST_AUTO_TEST_CASE (OverlapsResolver_two_overlaps) {
    using namespace npge;
    /*
    Input:
        Block 1:
            seq0: ---xxxx----
            seq1: ---xxxx----
            seq2: ---xxxx----
        Block 2
            seq1: -----xxxx--
            seq2: -----xxxx--
            seq3: -----xxxx--
    Output of OverlapsResolver:
        Block 1:
            seq0: ---xx------
            seq1: ---xx------
            seq2: ---xx------
        Block 2
            seq1: -------xx--
            seq2: -------xx--
            seq3: -------xx--
        Block 3:
            seq0: -----xx----
            seq1: -----xx----
            seq2: -----xx----
            seq3: -----xx----
    */
    SequencePtr s0 = boost::make_shared<InMemorySequence>("ctgc|ACAG|gacgt");
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("ctgcac|AGGA|cgt");
    Fragment* f10 = new Fragment(s0, 4, 7, -1);
    Fragment* f11 = new Fragment(s1, 4, 7, -1);
    Fragment* f12 = new Fragment(s2, 4, 7, -1);
    Block* b1 = new Block();
    b1->insert(f10);
    b1->insert(f11);
    b1->insert(f12);
    Fragment* f21 = new Fragment(s1, 6, 9);
    Fragment* f22 = new Fragment(s2, 6, 9);
    Fragment* f23 = new Fragment(s3, 6, 9);
    Block* b2 = new Block();
    b2->insert(f21);
    b2->insert(f22);
    b2->insert(f23);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    OverlapsResolver resolver;
    resolver.apply(block_set);
    BOOST_REQUIRE(block_set->size() == 3);
    int b[5] = {0, 0, 0, 0, 0};
    BOOST_FOREACH (Block* block, *block_set) {
        BOOST_REQUIRE(block->size() <= 4);
        b[block->size()] += 1;
        if (block->size() == 3) {
            BOOST_CHECK(block->front()->str() == "AC" ||
                        block->front()->str() == "GT" ||
                        block->front()->str() == "GA" ||
                        block->front()->str() == "TC");
        } else if (block->size() == 4) {
            BOOST_CHECK(block->front()->str() == "AG" ||
                        block->front()->str() == "CT");
        } else {
            BOOST_ERROR("Bad block size");
        }
    }
    BOOST_CHECK(!b[0] && !b[1] && !b[2] && b[3] == 2 && b[4] == 1);
}

BOOST_AUTO_TEST_CASE (OverlapsResolver_internal_subfragment) {
    using namespace npge;
    /*
    Input:
        Block 1:
            seq0: -----xx----
            seq1: -----xx----
            seq2: -----xx----
        Block 2
            seq1: ---xxxxxx--
            seq2: ---xxxxxx--
            seq3: ---xxxxxx--
    Output of OverlapsResolver:
        Block 1:
            seq1: ---xx------
            seq2: ---xx------
            seq3: ---xx------
        Block 2
            seq1: -------xx--
            seq2: -------xx--
            seq3: -------xx--
        Block 3:
            seq0: -----xx----
            seq1: -----xx----
            seq2: -----xx----
            seq3: -----xx----
    */
    SequencePtr s0 = boost::make_shared<InMemorySequence>("ctgcac|AG|gacgt");
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("ctgc|ACAGGA|cgt");
    Fragment* f10 = new Fragment(s0, 6, 7, 1);
    Fragment* f11 = new Fragment(s1, 6, 7, 1);
    Fragment* f12 = new Fragment(s2, 6, 7, 1);
    Block* b1 = new Block();
    b1->insert(f10);
    b1->insert(f11);
    b1->insert(f12);
    Fragment* f21 = new Fragment(s1, 4, 9);
    Fragment* f22 = new Fragment(s2, 4, 9);
    Fragment* f23 = new Fragment(s3, 4, 9);
    Block* b2 = new Block();
    b2->insert(f21);
    b2->insert(f22);
    b2->insert(f23);
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    OverlapsResolver resolver;
    resolver.apply(block_set);
    BOOST_REQUIRE(block_set->size() == 3);
    int b[5] = {0, 0, 0, 0, 0};
    BOOST_FOREACH (Block* block, *block_set) {
        BOOST_REQUIRE(block->size() <= 4);
        b[block->size()] += 1;
        if (block->size() == 3) {
            BOOST_CHECK(block->front()->str() == "AC" ||
                        block->front()->str() == "GT" ||
                        block->front()->str() == "GA" ||
                        block->front()->str() == "TC");
            BOOST_CHECK(!block->front()->next() ||
                        block->front()->next()->block()->size() == 4);
            BOOST_CHECK(!block->front()->prev() ||
                        block->front()->prev()->block()->size() == 4);
        } else if (block->size() == 4) {
            BOOST_CHECK(block->front()->str() == "AG" ||
                        block->front()->str() == "CT");
            BOOST_CHECK(!block->front()->next() ||
                        block->front()->next()->block()->size() == 3);
            BOOST_CHECK(!block->front()->prev() ||
                        block->front()->prev()->block()->size() == 3);
        } else {
            BOOST_ERROR("Bad block size");
        }
    }
    BOOST_CHECK(!b[0] && !b[1] && !b[2] && b[3] == 2 && b[4] == 1);
}

BOOST_AUTO_TEST_CASE (OverlapsResolver_multioverlaps) {
    using namespace npge;
    SequencePtr s[10];
    for (int j = 0; j < 10; j++) {
        s[j] = boost::make_shared<InMemorySequence>("CTGCACAGGACGTTGCACGGACGT");
    }
    BlockSetPtr block_set = new_bs();
    for (int i = 0; i < 10; i++) {
        Block* b = new Block();
        for (int j = 0; j < 10; j++) {
            b->insert(new Fragment(s[j], i, i + 10));
        }
        block_set->insert(b);
    }
    Connector connector;
    connector.apply(block_set);
    OverlapsResolver resolver;
    resolver.apply(block_set);
}

BOOST_AUTO_TEST_CASE (OverlapsResolver_TACG) {
    using namespace npge;
    std::string TACG;
    while (TACG.length() < 518) {
        TACG += "TACG";
    }
    SequencePtr s1 = boost::make_shared<InMemorySequence>(TACG);
    OverlapsResolver resolver;
    resolver.block_set()->add_sequence(s1);
    BlockSet& bs = *resolver.block_set();
    Block* block = new Block;
    bs.insert(block);
    block->insert(new Fragment(s1, 7, 130, 1));
    block->insert(new Fragment(s1, 7, 130, -1));
    block->insert(new Fragment(s1, 135, 258, 1));
    block->insert(new Fragment(s1, 135, 258, -1));
    block->insert(new Fragment(s1, 263, 386, 1));
    block->insert(new Fragment(s1, 263, 386, -1));
    block->insert(new Fragment(s1, 391, 514, 1));
    block->insert(new Fragment(s1, 391, 514, -1));
    Connector connector;
    connector.apply(resolver.block_set());
    BOOST_CHECK(resolver.overlaps());
}

