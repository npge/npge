/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ANNOTATION_HPP_
#define NPGE_ANNOTATION_HPP_

#include <string>

namespace npge {

/** Returns if a line is ID-line in annotation file */
bool is_id(const std::string& line);

}

#endif
