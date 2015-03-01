/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/test/unit_test.hpp>

#include "CutGaps.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"

using namespace npge;

const char* CUT_GAPS_INPUT =
    ">a_0_10 block=a\n"
    "AAAA---AAAA---AAA\n"
    ">a_0_4 block=a\n"
    "------AAAAA------\n"
    ">a_0_6 block=a\n"
    "-----AAAAAAA-----\n"
    ">a_0_8 block=a\n"
    "AAA----AAAAAA----\n\n";

const char* CUT_GAPS_OUTPUT =
    ">a_0_4 block=a\n"
    "AAAAA\n"
    ">a_1_5 block=a\n"
    "AAAAA\n"
    ">a_3_6 block=a\n"
    "-AAAA\n"
    ">a_4_7 block=a\n"
    "-AAAA\n\n";

BOOST_AUTO_TEST_CASE (CutGaps_main) {
    CutGaps cut_gaps;
    SequencePtr s = boost::make_shared<InMemorySequence>("AAAAAAAAAAAAAAAAAA");
    s->set_name("A");
    cut_gaps.block_set()->add_sequence(s);
    std::stringstream input(CUT_GAPS_INPUT);
    input >> *cut_gaps.block_set();
    cut_gaps.run();
    std::stringstream output;
    output << *cut_gaps.block_set();
    BOOST_CHECK(output.str() == CUT_GAPS_OUTPUT);
}

