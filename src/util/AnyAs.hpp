/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ANYAS_HPP_
#define BR_ANYAS_HPP_

#include <boost/any.hpp>

namespace bloomrepeats {

/** Class for boost::any + .as<T> method */
class AnyAs : public boost::any {
public:
    /** Constructor */
    AnyAs()
    { }

    /** Constructor */
    template<typename ValueType>
    AnyAs(const ValueType& value):
        boost::any(value)
    { }

    /** Constructor */
    AnyAs(const boost::any& value):
        boost::any(value)
    { }

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
};

}

#endif

