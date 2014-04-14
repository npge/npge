/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstring>
#include <stdexcept>
#include <sstream>
#include "throw_assert.hpp"

namespace boost {

#define SRC_PATTERN "src/"

const char* reduce_path(const char* file) {
    const char* subpath = strstr(file, SRC_PATTERN);
    if (subpath) {
        subpath += strlen(SRC_PATTERN);
        return subpath;
    } else {
        return file;
    }
}

void assertion_failed(char const* expr, char const* function,
                      char const* file, long line) {
    std::stringstream err;
    err << reduce_path(file) << ":" << line << ": " << function << ": ";
    err << "Assertation `" << expr << "' failed.";
    throw std::logic_error(err.str());
}

void assertion_failed_msg(char const* expr, char const* msg,
                          char const* function,
                          char const* file, long line) {
    std::stringstream err;
    err << reduce_path(file) << ":" << line << ": " << function << ": ";
    err << "Assertation `" << expr << "' failed." << std::endl;
    err << "Error message `" << msg << "'.";
    throw std::logic_error(err.str());
}

}

