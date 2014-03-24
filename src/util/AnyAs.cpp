/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <vector>
#include <typeinfo>

#include "AnyAs.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

bool any_equal(const AnyAs& a, const AnyAs& b) {
    BOOST_ASSERT(good_opt_type(a.type()));
    BOOST_ASSERT(good_opt_type(b.type()));
    if (a.type() != b.type()) {
        return false;
    }
    if (a.type() == typeid(int)) {
        return a.as<int>() == b.as<int>();
    } else if (a.type() == typeid(bool)) {
        return a.as<bool>() == b.as<bool>();
    } else if (a.type() == typeid(double)) {
        return a.as<double>() == b.as<double>();
    } else if (a.type() == typeid(std::string)) {
        return a.as<std::string>() == b.as<std::string>();
    } else if (a.type() == typeid(std::vector<std::string>)) {
        return a.as<std::vector<std::string> >() ==
               b.as<std::vector<std::string> >();
    }
    throw Exception("wrong type of any");
}

bool good_opt_type(const std::type_info& ti) {
    return ti == typeid(int) || ti == typeid(bool) || ti == typeid(double) ||
           ti == typeid(std::string) || ti == typeid(std::vector<std::string>);
}

}

