/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_READ_FILE_HPP_
#define NPGE_READ_FILE_HPP_

#include <iosfwd>
#include <string>

namespace npge {

/** Read stream to string */
std::string read_stream(std::istream& stream);

/** Return constents of file */
std::string read_file(const std::string& filename);

/** Read std::cin */
std::string read_stdin();

}

#endif

