/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_STEM_HPP_
#define BR_STEM_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Filter out blocks not represented in at least one of genomes.
Genome name is obtained by Sequence::genome().
*/
class Stem : public Processor {
protected:
    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

