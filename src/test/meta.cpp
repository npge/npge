/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/test/unit_test.hpp>

#include "Meta.hpp"
#include "Processor.hpp"

using namespace npge;

BOOST_AUTO_TEST_CASE (Meta_main) {
    Meta m;
    BOOST_FOREACH (std::string key, m.keys()) {
        BOOST_CHECK(m.has(key));
        BOOST_CHECK(m.get(key)->key() == key);
    }
}

BOOST_AUTO_TEST_CASE (Meta_options) {
    Meta m;
    m.set_opt("TEST", 33);
    SharedProcessor p = m.get("Processor");
    p->add_gopt("test", "Test", "TEST");
    BOOST_CHECK(p->opt_value("test").as<int>() == 33);
    m.set_opt("TEST", 34);
    BOOST_CHECK(p->opt_value("test").as<int>() == 34);
    p->set_opt_value("test", 35);
    BOOST_CHECK(p->opt_value("test").as<int>() == 35);
    p->set_opt_value("test", std::string("$TEST"));
    BOOST_CHECK(p->opt_value("test").as<int>() == 34);
    m.set_opt("FOO", 40);
    p->set_options("--test=$FOO");
    BOOST_CHECK(p->opt_value("test").as<int>() == 40);
}

