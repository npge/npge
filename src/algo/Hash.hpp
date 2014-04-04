/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_HASH_HPP_
#define BR_HASH_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Writes hash of blockset to std::cerr. FIXME
*/
class Hash : public Processor {
public:
    /** Constructor */
    Hash();

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

