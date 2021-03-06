/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ANYAS_HPP_
#define NPGE_ANYAS_HPP_

#include <iosfwd>
#include <boost/any.hpp>

namespace std {
struct type_info;
}

namespace npge {

/** Class for boost::any + .as<T> method */
class AnyAs : public boost::any {
public:
    /** Constructor */
    AnyAs() {
    }

    /** Constructor */
    template<typename ValueType>
    AnyAs(const ValueType& value):
        boost::any(value) {
    }

    /** Constructor */
    AnyAs(const boost::any& value):
        boost::any(value) {
    }

    /** Shortcut for any_cast */
    template<typename T>
    const T& as() const {
        return boost::any_cast<const T&>(*this);
    }

    /** Shortcut for any_cast */
    template<typename T>
    T& as() {
        return boost::any_cast<T&>(*this);
    }

    /** Convert value to string.
    \warning Option must be of good_opt_type.
    */
    std::string to_s() const;

    /** Convert value to Lua expression,
    \warning Option must be of good_opt_type.
    */
    std::string to_lua() const;

    /** Read value from string.
    \warning Option must be of good_opt_type.
    \warning Option must not be empty.
    \warning Throws exceptions on bad lexical cast.
    */
    void from_s(const std::string& value);

    /** Return string representation of type */
    std::string type_name() const;
};

/** Compare two any values.
Any's must be of one of fillowing types:
bool, int, Decimal, string, vector<string>.
*/
bool any_equal(const AnyAs& a, const AnyAs& b);

/** Return if type of the option is good.
Any's must be of one of fillowing types:
bool, int, Decimal, string, vector<string>.
*/
bool good_opt_type(const std::type_info& ti);

/** Streaming operator */
std::ostream& operator<<(std::ostream& o, const AnyAs& a);

}

#endif

