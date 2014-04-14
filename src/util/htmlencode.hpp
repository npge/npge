/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_HTMLENCODE_HPP_
#define BR_HTMLENCODE_HPP_

#include <string>

namespace bloomrepeats {

/** Return input string with escaped HTML special chars */
std::string htmlencode(const std::string& text);

}

#endif

