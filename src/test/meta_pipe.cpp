/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "meta_pipe.hpp"
#include "Pipe.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"

using namespace bloomrepeats;

BOOST_AUTO_TEST_CASE (MetaPipe_main) {
    BOOST_CHECK(create_pipe("pipe Empty{};")->key() == "Empty");
    BOOST_CHECK(create_pipe("pipe Empty {};")->key() == "Empty");
    BOOST_CHECK(create_pipe("pipe Empty {};")->name() == "Empty");
    BOOST_CHECK(create_pipe("pipe Empty { name \"123\";};")->key() == "Empty");
    BOOST_CHECK(create_pipe("pipe Empty { name \"123\";};")->name() == "123");
    BOOST_CHECK(create_pipe("pipe Empty { workers 100;};")->workers() == 100);
    ProcessorPtr with_timing = create_pipe("pipe Empty { timing true;};");
    BOOST_CHECK(with_timing->timing() == true);
    with_timing->set_timing(false); // not to write to std err
    BOOST_CHECK(create_pipe("pipe Empty { timing false;};")->timing() == 0);
    BOOST_CHECK(create_pipe("pipe E{no_options true;};")->no_options() == 1);
    BOOST_CHECK(create_pipe("pipe E{no_options false;};")->no_options() == 0);
    BOOST_CHECK(create_pipe("pipe Empty { max_loops 50;};")->max_loops() == 50);
}

BOOST_AUTO_TEST_CASE (MetaPipe_add) {
    ProcessorPtr filter = create_pipe("pipe F { add Filter target=target;};");
    filter->block_set()->insert(new Block);
    filter->run();
    BOOST_CHECK(filter->block_set()->empty());
    filter = create_pipe("pipe F { add Filter;};");
    filter->block_set()->insert(new Block);
    filter->run();
    BOOST_CHECK(filter->block_set()->empty());
}

