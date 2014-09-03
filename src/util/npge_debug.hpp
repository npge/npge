/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_DEBUG_HPP_
#define NPGE_DEBUG_HPP_

namespace npge {

/** Enable/disable Debug mode for asserts.
In debug mode, assertation fail causes program crash
to preserve stack state.
*/
void set_npge_debug(bool debug);

}

#endif

