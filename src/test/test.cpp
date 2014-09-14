/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "Meta.hpp"
#include "throw_assert.hpp"
#include "name_to_stream.hpp"

namespace npge {

struct ArgsFixture {
    ArgsFixture() {
        using namespace boost::unit_test;
        set_app_path(framework::master_test_suite().argv[0]);
    }
};

BOOST_GLOBAL_FIXTURE(ArgsFixture);

BOOST_GLOBAL_FIXTURE(Meta);

}

