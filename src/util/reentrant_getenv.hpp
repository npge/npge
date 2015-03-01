/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REENTRANT_GETENV_HPP_
#define NPGE_REENTRANT_GETENV_HPP_

#include <string>

namespace npge {

/** Reentrant getenv.
System getenv is called under mutex.
This can still crash if other code uses dangerous getenv
simultaneously.
*/
std::string reentrant_getenv(const std::string& name);

}

#endif

