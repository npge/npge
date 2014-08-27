/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_TSS_META_HPP_
#define NPGE_TSS_META_HPP_

#include "global.hpp"
#include "AnyAs.hpp"

namespace npge {

/** Set thread-local instance of Meta.
Restores previous value in destructor.

It is used in constructors of processors,
until Processor.set_meta() is called.
*/
struct TssMetaHolder {
    /** Constructor */
    TssMetaHolder(Meta* meta);

    /** Destructor */
    ~TssMetaHolder();

    Meta* prev_;
};

/** Return default instance of Meta for current thread */
Meta* tss_meta();

}

#endif

