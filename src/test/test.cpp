/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/shared_ptr.hpp>

#include "Meta.hpp"
#include "throw_assert.hpp"
#include "name_to_stream.hpp"

namespace npge {

struct TestFixture {
    TestFixture() {
        using namespace boost::unit_test;
        set_app_path(framework::master_test_suite().argv[0]);
        meta_.reset(new Meta);
    }

    boost::shared_ptr<Meta> meta_;
};

BOOST_GLOBAL_FIXTURE(TestFixture);

}

