/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/test/unit_test.hpp>

#include "MoveGaps.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"
#include "Decimal.hpp"

using namespace npge;

const char* MOVE_GAPS_INPUT =
    ">a_0_4 block=a\n"
    "A------AAAA------\n"
    ">a_0_6 block=a\n"
    "AA-----AAAA-----A\n"
    ">a_0_8 block=a\n"
    "AAA----AAAA----AA\n"
    ">a_0_10 block=a\n"
    "AAAA---AAAA---AAA\n"
    ">a_0_3 block=a\n"
    "AAAA-------------\n\n";

const char* MOVE_GAPS_OUTPUT =
    ">a_0_10 block=a\n"
    "AAAA---AAAA---AAA\n"
    ">a_0_3 block=a\n"
    "AAAA-------------\n"
    ">a_0_4 block=a\n"
    "------AAAAA------\n"
    ">a_0_6 block=a\n"
    "-----AAAAAAA-----\n"
    ">a_0_8 block=a\n"
    "AAA----AAAAAA----\n\n";

BOOST_AUTO_TEST_CASE (MoveGaps_main) {
    MoveGaps move_gaps;
    move_gaps.set_opt_value("max-tail", 3);
    move_gaps.set_opt_value("max-tail-to-gap", D(0.5));
    SequencePtr s = boost::make_shared<InMemorySequence>("AAAAAAAAAAAAAAAAAA");
    s->set_name("A");
    move_gaps.block_set()->add_sequence(s);
    std::stringstream input(MOVE_GAPS_INPUT);
    input >> *move_gaps.block_set();
    move_gaps.run();
    std::stringstream output;
    output << *move_gaps.block_set();
    BOOST_CHECK(output.str() == MOVE_GAPS_OUTPUT);
}

