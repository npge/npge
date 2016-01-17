/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_VERSION_HPP_
#define NPGE_VERSION_HPP_

#include <lua.hpp>

namespace npge {

/** Adds global table npge and fields VERSION, ARCH, COMMIT */
void init_npge_version(lua_State* L);

}

#endif

