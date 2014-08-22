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
#include "Decimal.hpp"
#include "global.hpp"
#include "cast.hpp"

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
    } else if (a.type() == typeid(Decimal)) {
        return a.as<Decimal>() == b.as<Decimal>();
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
    } else if (type() == typeid(Decimal)) {
        return as<Decimal>().to_s();
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
    } else if (type() == typeid(Decimal)) {
        *this = Decimal(value);
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
    return ti == typeid(int) || ti == typeid(bool) ||
           ti == typeid(Decimal) ||
           ti == typeid(std::string) || ti == typeid(Strings);
}

std::ostream& operator<<(std::ostream& o, const AnyAs& a) {
    o << a.to_s();
    return o;
}

}

