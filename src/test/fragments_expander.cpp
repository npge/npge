/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "FragmentsExpander.hpp"
#include "Connector.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

BOOST_AUTO_TEST_CASE (FragmentsExpander_expand_basic) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b = new Block();
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    b->insert(f1);
    b->insert(f2);
    FragmentsExpander fragments_expander;
    fragments_expander.expand(b);
    BOOST_CHECK(f1->min_pos() == 0);
    BOOST_CHECK(f1->max_pos() == s1->size() - 1);
    BOOST_CHECK(f2->min_pos() == 0);
    BOOST_CHECK(f2->max_pos() == s2->size() - 1);
    delete b;
}

BOOST_AUTO_TEST_CASE (FragmentsExpander_expand_3) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAaatcagatcg");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGGTCCGAgcggacggcc");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TGGTCCGAgcggacggcc");
    Block* b = new Block();
    Fragment* f1 = new Fragment(s1, 1, 2);
    Fragment* f2 = new Fragment(s2, 1, 2);
    Fragment* f3 = new Fragment(s3, 1, 2);
    b->insert(f1);
    b->insert(f2);
    b->insert(f3);
    FragmentsExpander fragments_expander;
    fragments_expander.set_opt_value("max-errors", 1);
    fragments_expander.set_opt_value("gap-range", 0);
    fragments_expander.expand(b);
    BOOST_CHECK(f1->min_pos() == 0);
    BOOST_CHECK(f1->max_pos() == 7);
    BOOST_CHECK(f2->min_pos() == 0);
    BOOST_CHECK(f2->max_pos() == 7);
    BOOST_CHECK(f3->min_pos() == 0);
    BOOST_CHECK(f3->max_pos() == 7);
    delete b;
}

BOOST_AUTO_TEST_CASE (FragmentsExpander_expand_two_blocks) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccGAcggccgcgga");
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
    FragmentsExpander fragments_expander;
    fragments_expander.expand(b1); // overlap = 0
    BOOST_CHECK(f11->min_pos() == 0);
    BOOST_CHECK(f11->max_pos() == 5);
    BOOST_CHECK(f12->min_pos() == 0);
    BOOST_CHECK(f12->max_pos() == 5);
    fragments_expander.expand(b2); // overlap = 0
    BOOST_CHECK(f21->min_pos() == 11);
    BOOST_CHECK(f21->max_pos() == 17);
    BOOST_CHECK(f22->min_pos() == 6);
    BOOST_CHECK(f22->max_pos() == 12);
    fragments_expander.set_opt_value("max-errors", 0);
    fragments_expander.set_max_overlap(1);
    fragments_expander.expand(b1);
    BOOST_CHECK(f11->min_pos() == 0);
    BOOST_CHECK(f11->max_pos() == 6);
    BOOST_CHECK(f12->min_pos() == 0);
    BOOST_CHECK(f12->max_pos() == 6);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (FragmentsExpander_expand_overlap) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGTCCGAGCGGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGTCCGAGCGGAcggcc");
    Block* b1 = new Block();
    Fragment* f11 = new Fragment(s1, 1, 5);
    Fragment* f12 = new Fragment(s2, 1, 5);
    b1->insert(f11);
    b1->insert(f12);
    Block* b2 = new Block();
    Fragment* f21 = new Fragment(s1, 3, 12);
    Fragment* f22 = new Fragment(s2, 3, 12);
    b2->insert(f21);
    b2->insert(f22);
    Fragment::connect(f11, f21);
    Fragment::connect(f12, f22);
    FragmentsExpander fragments_expander;
    fragments_expander.expand(b1); // overlap = 0
    BOOST_CHECK(f11->min_pos() == 0);
    BOOST_CHECK(f11->max_pos() == 5);
    BOOST_CHECK(f12->min_pos() == 0);
    BOOST_CHECK(f12->max_pos() == 5);
    fragments_expander.expand(b2); // overlap = 0
    BOOST_CHECK(f21->min_pos() == 3);
    BOOST_CHECK(f21->max_pos() == 17);
    BOOST_CHECK(f22->min_pos() == 3);
    BOOST_CHECK(f22->max_pos() == 17);
    delete b1;
    delete b2;
}

BOOST_AUTO_TEST_CASE (FragmentsExpander_block_set) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcgGAcggcc");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccGAcggccgcgga");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("tGGtccgacggccgcgga");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    b1->insert(new Fragment(s3, 1, 2));
    Block* b2 = new Block();
    b2->insert(new Fragment(s1, 11, 12));
    b2->insert(new Fragment(s2, 6, 7));
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Connector connector;
    connector.apply(block_set);
    FragmentsExpander fragments_expander;
    fragments_expander.apply(block_set);
    BOOST_CHECK(b1->front()->length() == 6);
    BOOST_CHECK(b1->front()->min_pos() == 0);
    BOOST_CHECK(b1->front()->str() == "TGGTCC");
    BOOST_CHECK(b2->front()->length() == 7);
    BOOST_CHECK(b2->front()->str() == "GACGGCC");
}

