/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "htmlencode.hpp"

namespace bloomrepeats {

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

