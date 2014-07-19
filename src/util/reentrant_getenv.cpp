/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include "boost-xtime.hpp"
#include <boost/thread/mutex.hpp>

namespace npge {

boost::mutex reentrant_getenv_mutex_;

std::string reentrant_getenv(const std::string& name) {
    boost::mutex::scoped_lock lock(reentrant_getenv_mutex_);
    const char* value = getenv(name.c_str());
    if (value) {
        return value;
    } else {
        return std::string();
    }
}

}

