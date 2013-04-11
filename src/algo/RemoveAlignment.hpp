/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_REMOVE_ALIGNMENT_HPP_
#define BR_REMOVE_ALIGNMENT_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Remove alignment rows of all fragments */
class RemoveAlignment : public Processor {
protected:
    bool run_impl() const;
};

}

#endif

