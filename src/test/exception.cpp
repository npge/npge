/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Exception.hpp"

BOOST_AUTO_TEST_CASE (Exception_main) {
    using namespace bloomrepeats;
    std::string error_message;
    try {
        throw Exception("Error message");
    } catch (Exception& e) {
        error_message = e.what();
    }
    BOOST_REQUIRE(error_message == "Error message");
}

