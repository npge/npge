/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/test/unit_test.hpp>

#include "Meta.hpp"
#include "Processor.hpp"

using namespace bloomrepeats;

BOOST_AUTO_TEST_CASE (Meta_main) {
    Meta m;
    BOOST_FOREACH (std::string key, m.keys()) {
        BOOST_CHECK(m.has(key));
        BOOST_CHECK(m.get(key)->key() == key);
    }
}

