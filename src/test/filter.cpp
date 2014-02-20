/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/scoped_ptr.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "Filter.hpp"
#include "SizeLimits.hpp"

BOOST_AUTO_TEST_CASE (Filter_good_block) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("AAA");
    Filter filter;
    allow_everything(&filter);
    filter.set_options("--min-fragment=100");
    filter.set_options("--min-block=2");
    boost::scoped_ptr<Block> block(new Block);
    block->insert(new Fragment(s1, 0, 0));
    block->insert(new Fragment(s1, 1, 1));
    BOOST_CHECK(!filter.is_good_block(block.get()));
    filter.set_options("--min-fragment=1");
    BOOST_CHECK(filter.is_good_block(block.get()));
}

