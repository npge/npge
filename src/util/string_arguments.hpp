/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_STRING_ARGUMENTS_HPP_
#define BR_STRING_ARGUMENTS_HPP_

#include <string>
#include <vector>

namespace bloomrepeats {

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
    StringToArgv(const std::vector<std::string>& arguments,
                 const char* dummy_app = "dummy_app");

    /** Constructor.
    \param dummy_app App name. It is added as first argument, if not 0.
    */
    StringToArgv(const char* dummy_app = "dummy_app");

    ~StringToArgv();

    /** Append argument.
    This method invalidates argc() and argv().
    */
    void add_argument(const std::string& argument);

    /** Return argc */
    int argc() const;

    /** Return argv */
    char** argv() const;

private:
    mutable std::vector<char*> argv_;
};

}

#endif

