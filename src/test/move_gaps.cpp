/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/test/unit_test.hpp>

#include "MoveGaps.hpp"
#include "BlockSet.hpp"
#include "Sequence.hpp"

using namespace bloomrepeats;

const char* MOVE_GAPS_INPUT =
    ">a_0_4 block=a\n"
    "a------aaaa------\n"
    ">a_0_6 block=a\n"
    "aa-----aaaa-----a\n"
    ">a_0_8 block=a\n"
    "aaa----aaaa----aa\n"
    ">a_0_10 block=a\n"
    "aaaa---aaaa---aaa\n"
    ">a_0_3 block=a\n"
    "aaaa-------------\n\n";

const char* MOVE_GAPS_OUTPUT =
    ">a_0_10 block=a\n"
    "aaaa---aaaa---aaa\n"
    ">a_0_3 block=a\n"
    "aaaa-------------\n"
    ">a_0_4 block=a\n"
    "------aaaaa------\n"
    ">a_0_6 block=a\n"
    "-----aaaaaaa-----\n"
    ">a_0_8 block=a\n"
    "aaa----aaaaaa----\n\n";

BOOST_AUTO_TEST_CASE (MoveGaps_main) {
    MoveGaps move_gaps(/* max-tail */ 3, /* max-tail-to-gap */ 0.5);
    SequencePtr s = boost::make_shared<InMemorySequence>("aaaaaaaaaaaaaaaaaa");
    s->set_name("a");
    move_gaps.block_set()->add_sequence(s);
    std::stringstream input(MOVE_GAPS_INPUT);
    input >> *move_gaps.block_set();
    move_gaps.run();
    std::stringstream output;
    output << *move_gaps.block_set();
    BOOST_CHECK(output.str() == MOVE_GAPS_OUTPUT);
}

