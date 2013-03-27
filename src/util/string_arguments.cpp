/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstring>
#include <string>
#include <vector>
#include <boost/foreach.hpp>

#include "string_arguments.hpp"

namespace bloomrepeats {

StringToArgv::StringToArgv(const std::vector<std::string>& arguments,
                           const char* dummy_app) {
    if (dummy_app) {
        add_argument(dummy_app);
    }
    BOOST_FOREACH (const std::string& argument, arguments) {
        add_argument(argument);
    }
}

StringToArgv::StringToArgv(const char* dummy_app) {
    if (dummy_app) {
        add_argument(dummy_app);
    }
}

StringToArgv::~StringToArgv() {
    BOOST_FOREACH (char* arg, argv_) {
        delete[] arg;
    }
    argv_.clear();
}

void StringToArgv::add_argument(const std::string& argument) {
    if (!argv_.empty()) {
        argv_.pop_back(); // 0-terminator
    }
    char* arg = new char[argument.size() + 1];
    std::strcpy(arg, argument.c_str());
    argv_.push_back(arg);
    argv_.push_back(0); // 0-terminator
}

int StringToArgv::argc() const {
    return argv_.size() - 1;
}

char** StringToArgv::argv() const {
    return &argv_[0];
}

}

