/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Exception.hpp"

namespace bloomrepeats {

Exception::Exception(const std::string& message):
    message_(message)
{ }

Exception::~Exception() throw()
{ }

const char* Exception::what() const throw() {
    return message_.c_str();
}

}

