/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OPTS_LIB_HPP_
#define NPGE_OPTS_LIB_HPP_

#include "global.hpp"

namespace npge {

/** Add global options to Meta */
void add_opts(Meta* meta);

/** Update opt value from environment variable.
Return whether option was successfully updated.
If new value can not be converted to type of target option,
return false.
*/
bool read_env(Meta* meta, const std::string& name);

/** Update all opts from environment variables */
void read_all_env(Meta* meta);

}

#endif

