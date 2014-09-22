/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <typeinfo>

#include "luabind-format-signature.hpp"
#include <luabind/luabind.hpp>
#include <luabind/operator.hpp>
#include <luabind/iterator_policy.hpp>
#include <luabind/object.hpp>

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
#include "name_to_stream.hpp"
#include "throw_assert.hpp"

namespace luabind {

typedef default_converter<npge::Strings> dcS;

int dcS::compute_score(lua_State* L, int index) {
    return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
}

static npge::Strings dcS_from(lua_State* L, int index) {
    npge::Strings result;
    lua_pushnil(L); // first key
    if (index < 0) {
        index -= 1;
    }
    ASSERT_EQ(lua_type(L, index), LUA_TTABLE);
    while (lua_next(L, index) != 0) {
        // uses 'key' (at index -2)
        // and 'value' (at index -1)
        lua_getglobal(L, "tostring");
        lua_pushvalue(L, -2); // value
        lua_call(L, 1, 1);
        result.push_back(lua_tostring(L, -1));
        lua_pop(L, 2); // remove two 'value's
        // keep 'key' for next iteration
    }
    return result;
}

npge::Strings dcS::from(lua_State* L, int index) {
    return dcS_from(L, index);
}

static void dcS_to(lua_State* L, const npge::Strings& strings) {
    lua_createtable(L, strings.size(), 0);
    for (int i = 0; i < strings.size(); i++) {
        const std::string& s = strings[i];
        lua_pushstring(L, s.c_str());
        lua_rawseti(L, -2, i + 1);
    }
}

void dcS::to(lua_State* L, const npge::Strings& strings) {
    dcS_to(L, strings);
}

typedef default_converter<npge::AnyAs> dcA;

int dcA::compute_score(lua_State* L, int index) {
    int t = lua_type(L, index);
    if (t == LUA_TNUMBER) {
        return 0;
    } else if (t == LUA_TBOOLEAN) {
        return 0;
    } else if (t == LUA_TSTRING) {
        return 0;
    } else if (t == LUA_TTABLE) {
        return 0;
    }
    return -1;
}

npge::AnyAs dcA::from(lua_State* L, int index) {
    int t = lua_type(L, index);
    if (t == LUA_TBOOLEAN) {
        return bool(lua_toboolean(L, index));
    } else if (t == LUA_TNUMBER) {
        return int(lua_tointeger(L, index));
    } else if (t == LUA_TUSERDATA) {
        luabind::object o(from_stack(L, index));
        return object_cast<npge::Decimal>(o);
    } else if (t == LUA_TSTRING) {
        return std::string(lua_tostring(L, index));
    } else if (t == LUA_TTABLE) {
        // Strings
        return dcS_from(L, index);
    }
    return npge::AnyAs();
}

void dcA::to(lua_State* L, const npge::AnyAs& a) {
    if (a.type() == typeid(bool)) {
        lua_pushboolean(L, a.as<bool>());
    } else if (a.type() == typeid(int)) {
        lua_pushinteger(L, a.as<int>());
    } else if (a.type() == typeid(npge::Decimal)) {
        luabind::object d(L, a.as<npge::Decimal>());
        d.push(L);
    } else if (a.type() == typeid(std::string)) {
        lua_pushstring(L, a.as<std::string>().c_str());
    } else if (a.type() == typeid(npge::Strings)) {
        dcS_to(L, a.as<npge::Strings>());
    }
}

#if LUABIND_INT64_MISSING

typedef default_converter<int64_t> dcI64;

int dcI64::compute_score(lua_State* L, int index) {
    int t = lua_type(L, index);
    if (t == LUA_TNUMBER) {
        return 0;
    }
    return -1;
}

int64_t dcI64::from(lua_State* L, int index) {
    return lua_tointeger(L, index);
}

void dcI64::to(lua_State* L, const int64_t& a) {
    lua_pushinteger(L, a);
}

#endif

}

namespace npge {

static int decimal_sub_point() {
    return Decimal::sub_point;
}

static int decimal_digits() {
    return Decimal::digits;
}

static int char_to_size_str(const std::string& c) {
    ASSERT_EQ(c.size(), 1);
    return char_to_size(c[0]);
}

static std::string size_to_char_str(int c) {
    return std::string(1, size_to_char(c));
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
    const char* start = text.c_str();
    if (ori == -1) {
        start += text.length() - 1;
    }
    return TO_S(make_hash(start, text.length(), ori));
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
           .def(constructor<const Decimal&>())
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

}

namespace npge {

void set_arg(lua_State* L, const Strings& a) {
    std::string args_lua = AnyAs(a).to_lua();
    std::string arg = "arg = " + args_lua;
    luaL_dostring(L, arg.c_str());
}

}

extern "C" int init_util_lua(lua_State* L) {
    using namespace luabind;
    using namespace npge;
    open(L);
    module(L) [
        register_decimal(),
        def("proportion", &proportion),
        def("char_to_size", &char_to_size_str),
        def("size_to_char", &size_to_char_str),
        def("complement_char", &complement_char),
        def("complement_str", &complement_str),
        def("complement_hash", &complement_hash_str),
        def("make_hash", &make_hash_str),
        def("make_hash", &make_hash_str1),
        def("reuse_hash", &reuse_hash_str),
        def("rand_name", &rand_name),
        def("make_seed", &make_seed),
        def("resolve_home_dir", &resolve_home_dir)
    ];
    return 0;
}

