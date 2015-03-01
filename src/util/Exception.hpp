/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_EXCEPTION_HPP_
#define NPGE_EXCEPTION_HPP_

#include <exception>
#include <string>

#include "global.hpp"

namespace npge {

/** Class for exceptions */
class Exception : public std::exception {
public:
    /** Constructor */
    Exception(const std::string& message);

    /** Destructor */
    ~Exception() throw();

    /** Error message */
    const char* what() const throw();

private:
    std::string message_;
};

}

#endif

