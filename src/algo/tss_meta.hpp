/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_TSS_META_HPP_
#define BR_TSS_META_HPP_

#include "global.hpp"
#include "AnyAs.hpp"

namespace bloomrepeats {

/** Return default instance of Meta for current thread.
\note This instance is not constant.
*/
Meta* tss_meta();

/** Return global option from thread local Meta */
AnyAs tss_go(const std::string& key, const AnyAs& dflt = 0);

}

#endif

