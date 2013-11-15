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

StringToArgv::StringToArgv(const StringToArgv& other) {
    argv_.push_back(0); // 0-terminator
    for (int i = 0; i < other.argc(); i++) {
        add_argument(other.argv_[i]);
    }
}

StringToArgv::StringToArgv(int argc, char** argv) {
    argv_.push_back(0); // 0-terminator
    for (int i = 0; i < argc; i++) {
        add_argument(argv[i]);
    }
}

StringToArgv::StringToArgv(const std::vector<std::string>& arguments,
                           const char* dummy_app) {
    argv_.push_back(0); // 0-terminator
    if (dummy_app) {
        add_argument(dummy_app);
    }
    BOOST_FOREACH (const std::string& argument, arguments) {
        add_argument(argument);
    }
}

StringToArgv::StringToArgv(const char* dummy_app) {
    argv_.push_back(0); // 0-terminator
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
    argv_.pop_back(); // 0-terminator
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

