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
#include "Joiner.hpp"
#include "Filter.hpp"
#include "Connector.hpp"

BOOST_AUTO_TEST_CASE (BlockSet_connect) {
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
    Connector connector;
    connector.apply(block_set);
    BOOST_CHECK(f1->next() == f2);
    BOOST_CHECK(f2->prev() == f1);
    BOOST_CHECK(f2->next() == f3);
    BOOST_CHECK(f3->prev() == f2);
}

BOOST_AUTO_TEST_CASE (BlockSet_filter) {
    using namespace npge;
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
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    block_set->insert(b3);
    Filter filter;
    filter.set_opt_value("min-fragment", 3);
    filter.set_opt_value("min-block", 1);
    filter.apply(block_set);
    BOOST_CHECK(block_set->size() == 1);
}

