/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "AnchorFinder.hpp"

BOOST_AUTO_TEST_CASE (AnchorFinder_main) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgGTCCGagCGGACggcc");
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_anchor_size(5);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1);
    if (block_set->size() == 1) {
        Fragment* f = block_set->front()->front();
        BOOST_CHECK(f->str() == "gtccg" || f->str() == "cggac");
    }
}

BOOST_AUTO_TEST_CASE (AnchorFinder_palindrome_elimination) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("atgcat");
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_anchor_size(6);
    anchor_finder.set_palindromes_elimination(true);
    BOOST_REQUIRE(anchor_finder.palindromes_elimination());
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 0);
    //
    block_set->clear();
    anchor_finder.set_palindromes_elimination(false);
    BOOST_REQUIRE(!anchor_finder.palindromes_elimination());
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_only_ori) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgGTCCGagCGGACggcc");
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_anchor_size(5);
    anchor_finder.set_only_ori(0);
    BOOST_REQUIRE(anchor_finder.only_ori() == 0);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1);
    //
    block_set->clear();
    anchor_finder.set_only_ori(1);
    BOOST_REQUIRE(anchor_finder.only_ori() == 1);
    anchor_finder.run();
    BOOST_CHECK(block_set->size() == 0);
    //
    block_set->clear();
    anchor_finder.set_only_ori(-1);
    BOOST_REQUIRE(anchor_finder.only_ori() == -1);
    anchor_finder.run();
    BOOST_CHECK(block_set->size() == 0);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_only_ori_3) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("gGTCCGaGTCCGtGTCCG");
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_anchor_size(5);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1 && block_set->front()->size() == 3);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_one_from_long_repeat) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("gaaagaaa");
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_anchor_size(3);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_several_sequences) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("gaaagaaa");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("gaaagaaa");
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->add_sequence(s1);
    block_set->add_sequence(s2);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_anchor_size(3);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1 && block_set->front()->size() == 4);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_two_workers) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("gaaagaaa");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("gaaagaaa");
    BlockSetPtr block_set = boost::make_shared<BlockSet>();
    block_set->add_sequence(s1);
    block_set->add_sequence(s2);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_anchor_size(3);
    anchor_finder.set_workers(2);
    anchor_finder.run();
    BOOST_WARN(block_set->size() >= 1 && block_set->front()->size() == 4);
}

