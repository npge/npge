/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Exception.hpp"
#include "simple_task.hpp"

BOOST_AUTO_TEST_CASE (Exception_main) {
    using namespace npge;
    std::string error_message;
    try {
        throw Exception("Error message");
    } catch (Exception& e) {
        error_message = e.what();
    }
    BOOST_CHECK(error_message == "Error message");
}

static void throwing_function() {
    using namespace npge;
    throw Exception("Test exception");
}

BOOST_AUTO_TEST_CASE (Exception_thread_group) {
    using namespace npge;
    Tasks tasks;
    tasks.push_back(throwing_function);
    tasks.push_back(throwing_function);
    std::string error_message;
    try {
        do_tasks(tasks_to_generator(tasks), /* workers */ 10);
    } catch (Exception& e) {
        error_message = e.what();
    }
    BOOST_CHECK(error_message == "Test exception");
}

