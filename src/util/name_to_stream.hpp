/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_NAME_TO_STREAM_HPP_
#define BR_NAME_TO_STREAM_HPP_

#include <iosfwd>
#include <string>
#include <boost/shared_ptr.hpp>

namespace bloomrepeats {

/** Return input stream for given filename.
If name is empty, return std::cin.

If name starts with ':', returns std::istringstream.

Otherwise returns std::ifstream.

Previous results are cached. To get them deleted/closed, call remove_istream().

Such names are accepted by FileWriter, FileReader and AbstractOutput.

This function is thread-safe.
*/
boost::shared_ptr<std::istream> name_to_istream(const std::string& name);

/** Return output stream for given filename.
If name is empty, return std::cout.

If name starts with ':', returns std::ostringstream.

Otherwise returns std::ofstream.

Previous results are cached. To get them deleted/closed, call remove_ostream().

Such names are accepted by FileWriter, FileReader and AbstractOutput.

This function is thread-safe.
*/
boost::shared_ptr<std::ostream> name_to_ostream(const std::string& name);

/** Remove stream, created/returned by name_to_istream().
Do not remove files being used. Remove after all manipulations with file.

This function is thread-safe.
*/
void remove_istream(const std::string& name);

/** Remove stream, created/returned by name_to_ostream().
Do not remove files being used. Remove after all manipulations with file.

This function is thread-safe.
*/
void remove_ostream(const std::string& name, bool remove_file);

}

#endif

