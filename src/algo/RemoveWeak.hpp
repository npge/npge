/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_REMOVE_WEAK_HPP_
#define BR_REMOVE_WEAK_HPP_

#include "Processor.hpp"

namespace npge {

/** Remove all weak blocks */
class RemoveWeak : public Processor {
public:
    /** Constructor */
    RemoveWeak();

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

