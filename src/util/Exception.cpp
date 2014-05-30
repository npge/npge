/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Exception.hpp"

namespace npge {

Exception::Exception(const std::string& message):
    message_(message) {
}

Exception::~Exception() throw() {
}

const char* Exception::what() const throw() {
    return message_.c_str();
}

}

