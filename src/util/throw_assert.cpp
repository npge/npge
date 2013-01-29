/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <stdexcept>
#include <sstream>
#include "throw_assert.hpp"

namespace boost {

void assertion_failed(char const* expr, char const* function,
                      char const* file, long line) {
    std::stringstream err;
    err << file << ":" << line << ": " << function << ": ";
    err << "Assertation `" << expr << "' failed.";
    throw std::logic_error(err.str());
}

void assertion_failed_msg(char const* expr, char const* msg,
                          char const* function,
                          char const* file, long line) {
    std::stringstream err;
    err << file << ":" << line << ": " << function << ": ";
    err << "Assertation `" << expr << "' failed." << std::endl;
    err << "Error message `" << msg << "'.";
    throw std::logic_error(err.str());
}

}

