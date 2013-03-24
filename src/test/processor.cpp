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

