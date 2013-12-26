/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include "Rest.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Filter.hpp"
#include "Connector.hpp"

BOOST_AUTO_TEST_CASE (Rest_main) {
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
    BlockSetPtr rest = new_bs();
    Rest r(block_set);
    r.apply(rest);
    BOOST_CHECK(rest->size() == 5);
    Filter filter;
    filter.set_min_block_size(1);
    filter.set_min_fragment_length(2);
    filter.apply(rest);
    BOOST_CHECK(rest->size() == 3);
    filter.set_min_fragment_length(6);
    filter.apply(rest);
    BOOST_CHECK(rest->size() == 2);
    filter.set_min_fragment_length(8);
    filter.apply(rest);
    BOOST_CHECK(rest->size() == 2);
    filter.set_min_fragment_length(9);
    filter.apply(rest);
    BOOST_CHECK(rest->size() == 1);
}

BOOST_AUTO_TEST_CASE (Rest_self) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("AAA");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 1, 1));
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    Connector connector;
    connector.apply(block_set);
    Rest r(block_set);
    r.apply(block_set);
    BOOST_CHECK(block_set->size() == 3);
}

