/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
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
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    Block* b = new Block();
    b->insert(f1);
    b->insert(f2);
    BlockSetPtr source = new_bs();
    source->insert(b);
    Union u(source);
    BlockSetPtr dest = new_bs();
    u.apply(dest);
    BOOST_REQUIRE(dest->size() == 1);
    BOOST_CHECK(dest->front()->size() == 2);
}

BOOST_AUTO_TEST_CASE (Union_clone_block) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Fragment* f2 = new Fragment(s1, 9, 13, -1);
    Block* block = new Block();
    block->insert(f2);
    Block* block_copy = Union::clone_block(block);
    BOOST_CHECK(block_copy->size() == 1);
    BOOST_CHECK(f2->block());
    BOOST_CHECK(!block_copy->has(f2));
    BOOST_CHECK(block_copy->front() != f2);
    BOOST_CHECK(*block_copy->front() == *f2);
    delete block_copy;
    delete block;
}

BOOST_AUTO_TEST_CASE (Union_clone_block_set) {
    using namespace npge;
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
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    BlockSetPtr block_set_copy =
        Union::clone_block_set(block_set);
    BOOST_CHECK(block_set_copy->size() == 3);
}

