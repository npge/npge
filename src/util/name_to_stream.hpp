/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_NAME_TO_STREAM_HPP_
#define NPGE_NAME_TO_STREAM_HPP_

#include <iosfwd>
#include <string>
#include <boost/shared_ptr.hpp>

#include "global.hpp"

namespace npge {

/** Return input stream for given filename.
If name is predefined (see set_istream()), returns corresponding stream.

If name starts with ':' or is empty, returns std::istringstream.

Otherwise returns std::ifstream.

Previous results are cached. To get them deleted/closed, call remove_istream().

Such names are accepted by FileWriter, FileReader and AbstractOutput.

This function is thread-safe.
*/
boost::shared_ptr<std::istream> name_to_istream(const std::string& name);

/** Associate input stream with given filename.
Predefined input streams (can be changed using this function):
 - '' => std::cin.
 - ':cin' => std::cin.
*/
void set_istream(const std::string& name, boost::shared_ptr<std::istream> s);

/** Remove association of filename with input stream */
void remove_istream(const std::string& name);

/** Return output stream for given filename.
If name is predefined (see set_ostream()), returns corresponding stream.

If name starts with ':' or is empty, returns std::ostringstream.

Otherwise returns std::ofstream.

Previous results are cached. To get them deleted/closed, call remove_ostream().

Such names are accepted by FileWriter, FileReader and AbstractOutput.

This function is thread-safe.
*/
boost::shared_ptr<std::ostream> name_to_ostream(const std::string& name);

/** Associate input stream with given filename.
Predefined input streams (can be changed using this function):
 - '' => std::cout.
 - ':cout' => std::cout.
 - ':cerr' => std::cerr.
 - ':null' => null output, as /dev/null.
*/
void set_ostream(const std::string& name, boost::shared_ptr<std::ostream> s);

/** Remove association of filename with output stream */
void remove_ostream(const std::string& name);

/** Sets same stringstream(contents) with set_isstream and set_osstream */
void set_sstream(const std::string& name, const std::string& contents = "");

/** Remove association of filename with input and output streams */
void remove_stream(const std::string& name);

/** Remove stream, created/returned by name_to_ostream().
Do not remove files being used. Remove after all manipulations with file.

This function is thread-safe.
*/
void remove_file(const std::string& name);

/** Return home directory.
Returns path to home directory in Windows and Unix,
if fails, returns dftl.
*/
std::string get_home_dir(const std::string& dftl = ".");

/** Return current directory */
std::string get_current_dir(const std::string& dftl = ".");

/** Set application path */
void set_app_path(const std::string& path);

/** Return application path */
std::string get_app_path();

/** Return application directory */
std::string get_app_dir();

/** Concatenate paths */
std::string cat_paths(const std::string& path1,
                      const std::string& path2);

/** system_complete */
std::string system_complete(const std::string& p);

/** Return if path exists */
bool file_exists(const std::string& path);

/** Return if path is directory */
bool is_dir(const std::string& path);

/** Return list of directory children */
Strings dir_children(const std::string& dir);

/** Copy file (overwrite) */
void copy_file(const std::string& src,
               const std::string& dst);

/** strip directory from filename */
std::string to_filename(const std::string& p);

/** Escape filename */
std::string escape_path(const std::string& str);

}

#endif

