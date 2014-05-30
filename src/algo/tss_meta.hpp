/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_TSS_META_HPP_
#define NPGE_TSS_META_HPP_

#include "global.hpp"
#include "AnyAs.hpp"

namespace npge {

/** Return default instance of Meta for current thread.
\note This instance is not constant.
*/
Meta* tss_meta();

/** Return global option from thread local Meta */
AnyAs tss_go(const std::string& key, const AnyAs& dflt = 0);

}

#endif

