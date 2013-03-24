/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <cstring>

#include "class_name.hpp"

namespace bloomrepeats {

std::string class_name(const char* type_info_name) {
    size_t length = std::strlen(type_info_name);
    const char* last = type_info_name + length - 1;
    if (*last == 'E') {
        last--;
    }
    const char* first = last;
    while (!isalpha(*first) && first > type_info_name) {
        first--;
    }
    while (isalpha(*first) && first > type_info_name) {
        first--;
    }
    if (first < last) {
        first++;
    }
    return std::string(first, last - first + 1);
}

}

