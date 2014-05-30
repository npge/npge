/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_HTMLENCODE_HPP_
#define BR_HTMLENCODE_HPP_

#include <string>

namespace npge {

/** Return input string with escaped HTML special chars */
std::string htmlencode(const std::string& text);

}

#endif

