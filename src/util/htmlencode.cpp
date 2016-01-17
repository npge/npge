/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "htmlencode.hpp"

namespace npge {

std::string htmlencode(const std::string& text) {
    std::string result;
    result.reserve(text.length());
    BOOST_FOREACH (char c, text) {
        if (c == '&') {
            result += "&amp;";
        } else if (c == '\'') {
            result += "&apos;";
        } else if (c == '"') {
            result += "&quot;";
        } else if (c == '<') {
            result += "&lt;";
        } else if (c == '>') {
            result += "&gt;";
        } else {
            result += c;
        }
    }
    return result;
}

}

