/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include "throw_assert.hpp"
#include "npge_debug.hpp"

namespace npge {

#define SRC_PATTERN "src/"

bool npge_debug_ = false;

const char* reduce_path(const char* file) {
    const char* subpath = strstr(file, SRC_PATTERN);
    if (subpath) {
        subpath += strlen(SRC_PATTERN);
        return subpath;
    } else {
        return file;
    }
}

void assertion_failed_msg(char const* expr, char const* msg,
                          char const* function,
                          char const* file, long line) {
    std::stringstream err;
    err << reduce_path(file) << ":" << line << ": " << function << ": ";
    err << "Assertation `" << expr << "' failed." << std::endl;
    err << "Error message `" << msg << "'.";
    if (npge_debug_) {
        std::cerr << err.str() << "\n";
        abort();
    } else {
        throw std::logic_error(err.str());
    }
}

void set_npge_debug(bool debug) {
    npge_debug_ = debug;
}

}

