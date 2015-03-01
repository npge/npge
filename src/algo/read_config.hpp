/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_READ_OPTS_HPP_
#define NPGE_READ_OPTS_HPP_

#include "global.hpp"

namespace npge {

/** Update opt value from environment variable.
Return whether option was successfully updated.
If new value can not be converted to type of target option,
throws.
*/
bool read_env(Meta* meta, const std::string& name);

/** Update all opts from environment variables */
void read_all_env(Meta* meta);

/** Read config files and environment variables.
Values from options CONFIG0..CONFIG9 are read
one after another. Empty values are skipped.
ENV means "read environment variables".
LOCAL_CONF means "read config file $LOCAL_CONF".
Other values are interpreted as file names.
Format of config file is Lua script.
*/
void read_config(Meta* meta);

}

#endif

