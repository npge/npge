/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_DOWNLOAD_FILE_HPP_
#define NPGE_DOWNLOAD_FILE_HPP_

#include <string>

namespace npge {

/** Download file using libCURL.
Return if success.
*/
bool download_file(const std::string& url,
                   const std::string& out_fname);

}

#endif
