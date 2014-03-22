/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MARK_NON_WEAK_HPP_
#define BR_MARK_NON_WEAK_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Mark all blocks non-weak.
\see Block::set_weak()

\note If Block share same fragments, then some of blocks are no remain weak.
*/
class MarkNonWeak : public Processor {
protected:
    void run_impl() const;
};

}

#endif

