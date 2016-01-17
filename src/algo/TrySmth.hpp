/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_TRY_SMTH_HPP
#define NPGE_TRY_SMTH_HPP

#include "Pipe.hpp"

namespace npge {

/** Try to do something, align (and filter), restore original if worse */
class TrySmth : public Pipe {
public:
    /** Constructor */
    TrySmth();

protected:
    const char* name_impl() const;
};

class AddingLoop;

/** Align and move overlapless from other to target.
Preserve non-overlapping parts of overlapping blocks,
align them again and so on.
*/
class AddingLoopBySize : public Processor {
public:
    /** Constructor */
    AddingLoopBySize();

protected:
    void run_impl() const;
    const char* name_impl() const;

private:
    AddingLoop* al_;
};

}

#endif

