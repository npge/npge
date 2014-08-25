/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>

#include "string_arguments.hpp"
#include "po.hpp"

namespace npge {

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

StringToArgv::StringToArgv(const Strings& arguments,
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

bool StringToArgv::remove_argument(const std::string& argument) {
    if (argv_.empty()) {
        return false;
    } else {
        int old_size = argv_.size();
        argv_.pop_back(); // 0-terminator
        std::vector<char*> new_argv;
        for (int i = 0; i < argv_.size(); i++) {
            if (argv_[i] == argument) {
                if (i + 1 < argv_.size() && argv_[i + 1]) {
                    if (argv_[i + 1][0] != '-') {
                        i++;
                    }
                }
            } else {
                new_argv.push_back(argv_[i]);
            }
        }
        argv_.swap(new_argv);
        argv_.push_back(0); // 0-terminator
        return argv_.size() < old_size;
    }
}

bool StringToArgv::has_argument(const std::string& argument) const {
    return has_arg(argc(), argv(), argument.c_str());
}

std::string StringToArgv::get_argument(const std::string& arg) const {
    bool eq = false;
    BOOST_FOREACH (const char* a, argv_) {
        if (a && a == arg) {
            eq = true;
        } else if (eq) {
            if (a && a[0] != '\0' && a[0] != '-') {
                return a;
            } else {
                return "";
            }
        }
    }
    return "";
}

int StringToArgv::argc() const {
    return argv_.size() - 1;
}

char** StringToArgv::argv() const {
    return &argv_[0];
}

std::string StringToArgv::to_s() const {
    std::string result;
    BOOST_FOREACH (const char* a, argv_) {
        if (a) {
            if (!result.empty()) {
                result += ' ';
            }
            // TODO escape
            result += a;
        }
    }
    return result;
}

Strings StringToArgv::to_strings() const {
    Strings result;
    BOOST_FOREACH (const char* a, argv_) {
        if (a) {
            result.push_back(a);
        }
    }
    return result;
}

}

