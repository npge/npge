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
#include "read_config.hpp"
#include "throw_assert.hpp"

namespace npge {

Meta meta;
TssMetaHolder tmh(&meta);

struct MetaReadConfig {
    MetaReadConfig() {
        ASSERT_EQ(tss_meta(), &meta);
        read_config(tss_meta());
    }
};

}

