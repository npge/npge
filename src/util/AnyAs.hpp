/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ANYAS_HPP_
#define BR_ANYAS_HPP_

#include <boost/any.hpp>

namespace std {
struct type_info;
}

namespace bloomrepeats {

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

    /** Read value from string.
    \warning Option must be of good_opt_type.
    \warning Option must not be empty.
    \warning Throws exceptions on bad lexical cast.
    */
    void from_s(const std::string& value);
};

/** Compare two any values.
Any's must be of one of fillowing types:
bool, int, double, string, vector<string>.
*/
bool any_equal(const AnyAs& a, const AnyAs& b);

/** Return if type of the option is good.
Any's must be of one of fillowing types:
bool, int, double, string, vector<string>.
*/
bool good_opt_type(const std::type_info& ti);

}

#endif

