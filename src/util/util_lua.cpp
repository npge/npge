/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <typeinfo>

#include <luabind/luabind.hpp>
#include <luabind/operator.hpp>
#include <luabind/iterator_policy.hpp>

#include "global.hpp"
#include "util_lua.hpp"
#include "proportion.hpp"
#include "cast.hpp"
#include "char_to_size.hpp"
#include "make_hash.hpp"
#include "complement.hpp"
#include "rand_name.hpp"
#include "Decimal.hpp"
#include "AnyAs.hpp"
#include "throw_assert.hpp"

namespace npge {

static int decimal_sub_point() {
    return Decimal::sub_point;
}

static int decimal_digits() {
    return Decimal::digits;
}

static bool anyas_as_bool(const AnyAs& a) {
    return a.as<bool>();
}

static int anyas_as_int(const AnyAs& a) {
    return a.as<int>();
}

static Decimal anyas_as_decimal(const AnyAs& a) {
    return a.as<Decimal>();
}

static std::string anyas_as_string(const AnyAs& a) {
    return a.as<std::string>();
}

static Strings anyas_as_strings(const AnyAs& a) {
    return a.as<Strings>();
}

static bool anyas_is_good(const AnyAs& a) {
    return good_opt_type(a.type());
}

static std::string anyas_type(const AnyAs& a) {
    if (a.type() == typeid(bool)) {
        return "bool";
    } else if (a.type() == typeid(int)) {
        return "int";
    } else if (a.type() == typeid(Decimal)) {
        return "Decimal";
    } else if (a.type() == typeid(std::string)) {
        return "std::string";
    } else if (a.type() == typeid(Strings)) {
        return "Strings";
    } else {
        return a.type().name();
    }
}

static std::string strings_at(const Strings& v, int index) {
    return v.at(index);
}

static const Strings& strings_iter(const Strings& v) {
    return v;
}

static std::string complement_char(const std::string c) {
    ASSERT_EQ(c.size(), 1);
    return std::string(1, complement(c[0]));
}

static std::string complement_str(const std::string& str) {
    std::string copy = str;
    complement(copy);
    return copy;
}

static std::string complement_hash_str(const std::string& hash0,
                                       int letters_number) {
    return TO_S(complement_hash(L_CAST<hash_t>(hash0),
                                letters_number));
}

static std::string make_hash_str(const std::string& text,
                                 int ori) {
    return TO_S(make_hash(text.c_str(), text.length(), ori));
}

static std::string make_hash_str1(const std::string& text) {
    return make_hash_str(text, 1);
}

static std::string reuse_hash_str(const std::string& old_hash0,
                                  size_t length,
                                  const std::string& rc,
                                  const std::string& ac,
                                  bool forward) {
    ASSERT_EQ(rc.size(), 1);
    ASSERT_EQ(ac.size(), 1);
    hash_t old_hash = L_CAST<hash_t>(old_hash0);
    hash_t new_hash = reuse_hash(old_hash, length,
                                 rc[0], ac[0], forward);
    return TO_S(new_hash);
}

luabind::scope register_decimal() {
    using namespace luabind;
    return class_<Decimal>("Decimal")
           .scope [
               def("sub_point", &decimal_sub_point),
               def("digits", &decimal_digits)
           ]
           .def(constructor<>())
           .def(constructor<int>())
           .def(constructor<const std::string&>())
           .def(const_self + const_self)
           .def(const_self - const_self)
           .def(-const_self)
           .def(const_self * const_self)
           .def(const_self / const_self)
           .def(const_self == const_self)
           .def(const_self < const_self)
           .def(const_self <= const_self)
           .def("to_d", &Decimal::to_d)
           .def("to_i", &Decimal::to_i)
           .def("fraction", &Decimal::fraction)
           .def("round", &Decimal::round)
           .def("to_s", &Decimal::to_s)
           .def(tostring(self))
          ;
}

luabind::scope register_anyas() {
    using namespace luabind;
    return class_<AnyAs>("AnyAs")
           .def(constructor<>())
           .def(constructor<bool>())
           .def(constructor<int>())
           .def(constructor<Decimal>())
           .def(constructor<std::string>())
           .def(constructor<Strings>())
           .def(tostring(self))
           .def("as_bool", &anyas_as_bool)
           .def("as_int", &anyas_as_int)
           .def("as_decimal", &anyas_as_decimal)
           .def("as_string", &anyas_as_string)
           .def("as_strings", &anyas_as_strings)
           .def("to_s", &AnyAs::to_s)
           .def("from_s", &AnyAs::from_s)
           .def("any_equal", &any_equal)
           .def("is_good", &anyas_is_good)
           .def("type", &anyas_type)
          ;
}

luabind::scope register_strings() {
    using namespace luabind;
    return class_<Strings>("Strings")
           .def(constructor<>())
           .def(constructor<int>())
           .def(constructor<int, std::string>())
           .def("empty", &Strings::empty)
           .def("clear", &Strings::clear)
           .def("size", &Strings::size)
           .def("resize", &Strings::resize)
           .def("push_back", &Strings::push_back)
           .def("at", &strings_at)
           .def("iter", &strings_iter, return_stl_iterator)
          ;
}

}

extern "C" int init_util_lua(lua_State* L) {
    using namespace luabind;
    using namespace npge;
    open(L);
    module(L) [
        register_decimal(),
        register_anyas(),
        register_strings(),
        def("proportion", &proportion),
        def("char_to_size", &char_to_size),
        def("size_to_char", &size_to_char),
        def("complement_char", &complement_char),
        def("complement_str", &complement_str),
        def("complement_hash", &complement_hash_str),
        def("make_hash", &make_hash_str),
        def("make_hash", &make_hash_str1),
        def("reuse_hash", &reuse_hash_str),
        def("rand_name", &rand_name),
        def("make_seed", &make_seed)
    ];
    return 0;
}

