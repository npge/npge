/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_STRING_ARGUMENTS_HPP_
#define NPGE_STRING_ARGUMENTS_HPP_

#include <vector>

#include "global.hpp"

namespace npge {

/** Convert std::string arguments to argc and argv */
class StringToArgv {
public:
    /** Constructor */
    StringToArgv(const StringToArgv& other);

    /** Constructor */
    StringToArgv(int argc, char** argv);

    /** Constructor.
    \param arguments Vector of arguemnts.
    \param dummy_app App name. It is added as first argument, if not 0.
    */
    StringToArgv(const Strings& arguments,
                 const char* dummy_app = "dummy_app");

    /** Constructor.
    \param dummy_app App name. It is added as first argument, if not 0.
    */
    StringToArgv(const char* dummy_app = "dummy_app");

    /** Destructor */
    ~StringToArgv();

    /** Append argument.
    This method invalidates argc() and argv().
    */
    void add_argument(const std::string& argument);

    /** Remove argument.
    Return if the argument was removed.
    */
    bool remove_argument(const std::string& argument);

    /** Return if the argument is in the argument list */
    bool has_argument(const std::string& argument) const;

    /** Return argument value or "" */
    std::string get_argument(const std::string& arg) const;

    /** Return argc */
    int argc() const;

    /** Return argv */
    char** argv() const;

    /** Join all arguments and return as string */
    std::string to_s() const;

    /** Join all arguments as list of strings*/
    Strings to_strings() const;

private:
    mutable std::vector<char*> argv_;
};

}

#endif

