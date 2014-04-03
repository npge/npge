/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_REMOVE_WEAK_HPP_
#define BR_REMOVE_WEAK_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Remove all weak blocks */
class RemoveWeak : public Processor {
protected:
    void run_impl() const;
};

}

#endif

