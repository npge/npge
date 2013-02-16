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
#include "StickBoundaries.hpp"

BOOST_AUTO_TEST_CASE (StickBoundaries_main) {
    using namespace bloomrepeats;
    std::string aaa(1000, 'a');
    SequencePtr s1 = boost::make_shared<InMemorySequence>(aaa);
    SequencePtr s2 = boost::make_shared<InMemorySequence>(aaa);
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    block_set->add_sequence(s2);
    Fragment* f11 = new Fragment(s1, 0, 999, -1);
    Fragment* f12 = new Fragment(s1, 10, 200);
    Fragment* f21 = new Fragment(s2, 0, 999, -1);
    Fragment* f21a = new Fragment(s2, 0, 979, -1);
    Fragment* f22 = new Fragment(s2, 100, 199);
    Fragment* f23 = new Fragment(s2, 100, 199, -1);
    Fragment* f24 = new Fragment(s2, 201, 299);
    Fragment* f25 = new Fragment(s2, 600, 749, -1);
    Fragment* f26 = new Fragment(s2, 600, 699);
    Block* b1 = new Block; // dummy
    b1->insert(f11);
    b1->insert(f12);
    b1->insert(f21);
    b1->insert(f21a);
    b1->insert(f22);
    b1->insert(f23);
    b1->insert(f24);
    b1->insert(f25);
    b1->insert(f26);
    block_set->insert(b1);
    StickBoundaries stick(/* min_distance */ 100);
    stick.set_block_set(block_set);
    stick.run();
    BOOST_CHECK(*f11 == Fragment(s1, 0, 999, -1));
    BOOST_CHECK(*f12 == Fragment(s1, 0, 200));
    BOOST_CHECK(*f21 == Fragment(s2, 0, 999, -1));
    BOOST_CHECK(*f21a == Fragment(s2, 0, 999, -1));
    BOOST_CHECK(*f22 == Fragment(s2, 100, 199));
    BOOST_CHECK(*f23 == Fragment(s2, 100, 199, -1));
    BOOST_CHECK(*f24 == Fragment(s2, 200, 299));
    BOOST_CHECK(*f25 == Fragment(s2, 600, 724, -1));
    BOOST_CHECK(*f26 == Fragment(s2, 600, 724));
}

