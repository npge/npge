/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "key_value.hpp"

namespace bloomrepeats {

std::string extract_value(const std::string& values, const std::string& key) {
    size_t key_pos = values.find(key + "=");
    if (key_pos == std::string::npos) {
        return "";
    }
    size_t value_start = key_pos + key.size() + 1; // 1 because of '='
    size_t space_pos = values.find(' ', value_start); // or npos
    return values.substr(value_start, space_pos - value_start);
}

}

