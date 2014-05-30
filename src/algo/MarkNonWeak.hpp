/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MARK_NON_WEAK_HPP_
#define BR_MARK_NON_WEAK_HPP_

#include "Processor.hpp"

namespace npge {

/** Mark all blocks non-weak.
\see Block::set_weak()

\note If Block share same fragments, then some of blocks
are now remain weak.
*/
class MarkNonWeak : public Processor {
public:
    /** Constructor */
    MarkNonWeak();

protected:
    void run_impl() const;
    const char* name_impl() const;
};

}

#endif

