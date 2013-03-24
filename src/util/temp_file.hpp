/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_TEMP_FILE_HPP_
#define BR_TEMP_FILE_HPP_

#include <string>

namespace bloomrepeats {

/** Return unique temp file name.
Boost.Filesystem's unique_path is used, if present.

Otherwise this function call tmpnam() from C++ Standard library,
checking the result to be writable.
If this fails 10 times, empty string is returned.
*/
std::string temp_file();

}

#endif

