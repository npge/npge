/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <string>
#include <vector>
#include <typeinfo>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>

#include "AnyAs.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "global.hpp"
#include "to_s.hpp"

namespace npge {

bool any_equal(const AnyAs& a, const AnyAs& b) {
    ASSERT_TRUE(good_opt_type(a.type()));
    ASSERT_TRUE(good_opt_type(b.type()));
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
    } else if (a.type() == typeid(Strings)) {
        return a.as<Strings>() ==
               b.as<Strings>();
    }
    throw Exception("wrong type of any");
}

std::string AnyAs::to_s() const {
    if (type() == typeid(bool)) {
        return TO_S(as<bool>());
    } else if (type() == typeid(int)) {
        return TO_S(as<int>());
    } else if (type() == typeid(double)) {
        return TO_S(as<double>());
    } else if (type() == typeid(std::string)) {
        return TO_S(as<std::string>());
    } else if (type() == typeid(Strings)) {
        using namespace boost::algorithm;
        return TO_S(join(as<Strings>(), " "));
    }
    std::string type_str = type().name();
    throw Exception("wrong type of any: " + type_str);
}

void AnyAs::from_s(const std::string& value) {
    ASSERT_FALSE(empty());
    if (type() == typeid(bool)) {
        *this = L_CAST<bool>(value);
        return;
    } else if (type() == typeid(int)) {
        *this = L_CAST<int>(value);
        return;
    } else if (type() == typeid(double)) {
        *this = L_CAST<double>(value);
        return;
    } else if (type() == typeid(std::string)) {
        *this = value;
        return;
    } else if (type() == typeid(Strings)) {
        using namespace boost::algorithm;
        Strings parts;
        split(parts, value, isspace, token_compress_on);
        *this = parts;
        return;
    }
    std::string type_str = type().name();
    throw Exception("wrong type of any: " + type_str);
}

bool good_opt_type(const std::type_info& ti) {
    return ti == typeid(int) || ti == typeid(bool) || ti == typeid(double) ||
           ti == typeid(std::string) || ti == typeid(Strings);
}

}

