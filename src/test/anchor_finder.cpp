/*
 * NPG-explorer, Nucleotide PanGenome explorer
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
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgGTCCGagCGGACggcc");
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 5);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1);
    if (block_set->size() == 1) {
        Fragment* f = block_set->front()->front();
        BOOST_CHECK(f->str() == "GTCCG" || f->str() == "CGGAC");
    }
}

BOOST_AUTO_TEST_CASE (AnchorFinder_n_negative) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgGTNCGagCGNACggcc");
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 5);
    anchor_finder.run();
    BOOST_CHECK(block_set->size() == 0);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_n_positive) {
    using namespace npge;
    std::string s = "GTNCGATAnnnGTNCGATA";
    SequencePtr s1 = boost::make_shared<InMemorySequence>(s);
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 5);
    anchor_finder.run();
    BOOST_CHECK(block_set->size() > 0);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_palindrome_elimination) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("ATGCAT");
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 6);
    anchor_finder.set_opt_value("no-palindromes", true);
    BOOST_REQUIRE(anchor_finder.opt_value("no-palindromes").as<bool>());
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 0);
    //
    block_set->clear();
    anchor_finder.set_opt_value("no-palindromes", false);
    BOOST_REQUIRE(!anchor_finder.opt_value("no-palindromes").as<bool>());
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_only_ori) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgGTCCGagCGGACggcc");
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 5);
    anchor_finder.set_opt_value("only-ori", 0);
    BOOST_REQUIRE(anchor_finder.opt_value("only-ori").as<int>() == 0);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1);
    //
    block_set->clear();
    anchor_finder.set_opt_value("only-ori", 1);
    BOOST_REQUIRE(anchor_finder.opt_value("only-ori").as<int>() == 1);
    anchor_finder.run();
    BOOST_CHECK(block_set->size() == 0);
    //
    block_set->clear();
    anchor_finder.set_opt_value("only-ori", -1);
    BOOST_REQUIRE(anchor_finder.opt_value("only-ori").as<int>() == -1);
    anchor_finder.run();
    BOOST_CHECK(block_set->size() == 0);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_only_ori_3) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("gGTCCGaGTCCGtGTCCG");
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 5);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1 && block_set->front()->size() == 3);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_one_from_long_repeat) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GAAAGAAA");
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 3);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_several_sequences) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GAAAGAAA");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("GAAAGAAA");
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    block_set->add_sequence(s2);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 3);
    anchor_finder.run();
    BOOST_WARN(block_set->size() == 1 && block_set->front()->size() == 4);
}

BOOST_AUTO_TEST_CASE (AnchorFinder_two_workers) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("GAAAGAAA");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("GAAAGAAA");
    BlockSetPtr block_set = new_bs();
    block_set->add_sequence(s1);
    block_set->add_sequence(s2);
    AnchorFinder anchor_finder;
    anchor_finder.set_block_set(block_set);
    anchor_finder.set_opt_value("anchor-size", 3);
    anchor_finder.set_workers(2);
    anchor_finder.run();
    BOOST_WARN(block_set->size() >= 1 && block_set->front()->size() == 4);
}

