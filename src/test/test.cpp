/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "Meta.hpp"
#include "tss_meta.hpp"
#include "throw_assert.hpp"

namespace npge {

struct TestMeta {
    Meta meta_;
    TssMetaHolder tmh_;

    TestMeta():
        tmh_(&meta_) {
    }
};

BOOST_GLOBAL_FIXTURE(TestMeta);

}

