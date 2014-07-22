/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "name_to_stream.hpp"
#include "read_file.hpp"

typedef boost::shared_ptr<std::ostream> OPtr;

BOOST_AUTO_TEST_CASE (name_to_stream_main) {
    using namespace npge;
    set_sstream(":o");
    OPtr o = name_to_ostream(":o");
    (*o) << "test";
    BOOST_CHECK(read_file(":o") == "test");
}

