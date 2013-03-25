/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/test/unit_test.hpp>

#include "CutGaps.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"

using namespace bloomrepeats;

const char* CUT_GAPS_INPUT =
    ">a_0_10 block=a\n"
    "aaaa---aaaa---aaa\n"
    ">a_0_4 block=a\n"
    "------aaaaa------\n"
    ">a_0_6 block=a\n"
    "-----aaaaaaa-----\n"
    ">a_0_8 block=a\n"
    "aaa----aaaaaa----\n\n";

const char* CUT_GAPS_OUTPUT =
    ">a_0_4 block=a\n"
    "aaaaa\n"
    ">a_1_5 block=a\n"
    "aaaaa\n"
    ">a_3_6 block=a\n"
    "-aaaa\n"
    ">a_4_7 block=a\n"
    "-aaaa\n\n";

BOOST_AUTO_TEST_CASE (CutGaps_main) {
    CutGaps cut_gaps;
    SequencePtr s = boost::make_shared<InMemorySequence>("aaaaaaaaaaaaaaaaaa");
    s->set_name("a");
    cut_gaps.block_set()->add_sequence(s);
    std::stringstream input(CUT_GAPS_INPUT);
    input >> *cut_gaps.block_set();
    cut_gaps.run();
    std::stringstream output;
    output << *cut_gaps.block_set();
    BOOST_CHECK(output.str() == CUT_GAPS_OUTPUT);
}

