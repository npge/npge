/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Processor.hpp"
#include "Filter.hpp"
#include "OverlapsResolver2.hpp"

using namespace bloomrepeats;

BOOST_AUTO_TEST_CASE (processor_name_main) {
    ProcessorPtr p(new Filter);
    BOOST_CHECK(processor_name(p) == "Filter");
    ProcessorPtr p1(new OverlapsResolver2);
    BOOST_CHECK(processor_name(p1) == "OverlapsResolver2");
}

BOOST_AUTO_TEST_CASE (processor_set_options) {
    Processor p1;
    BOOST_CHECK(!p1.no_options());
    BOOST_CHECK(!p1.timing());
    //
    Processor p2;
    p2.set_options("no_options");
    BOOST_CHECK(p2.no_options());
    //
    Processor p3;
    p3.set_options("--timing");
    BOOST_CHECK(p3.timing());
    p3.set_timing(false); // not to write to std err
}

