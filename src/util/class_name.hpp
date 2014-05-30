/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CLASS_NAME_HPP_
#define BR_CLASS_NAME_HPP_

#include <string>

namespace npge {

/** Return class name based on type_info::name().
Remove terminal 'E' if there.
Beginning from the end of string, find first letter,
then first non-letter (i).
Return slice [i + 1, end).
*/
std::string class_name(const char* type_info_name);

}

#endif

