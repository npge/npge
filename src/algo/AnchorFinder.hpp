/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ANCHOR_FINDER_HPP_
#define NPGE_ANCHOR_FINDER_HPP_

#include "global.hpp"
#include "Processor.hpp"

namespace npge {

/** Finder of short anchors.

For large repeats one short part is selected and
returned as a anchor (unless --anchor-similar=false).

AnchorFinder memorizes hashes of previous run()'s
and skips them from output.

\note Using >= 2 workers may (very unlikely) cause races,
    since bloom filter is not protected by a mutex.
    Such a races may cause some anchors not to be found.
\note The smallest piece of work, passed to a worker,
    is one sequence.
    So it is useless to set workers > sequences.

*/
class AnchorFinder : public Processor {
public:
    /** Default constructor */
    AnchorFinder();

    /** Destructor */
    ~AnchorFinder();

protected:
    /** Find anchors in added sequence.
    Each found anchor is inserted into block_set().
    */
    void run_impl() const;

    const char* name_impl() const;

private:
    class Impl;
    Impl* impl_;
};

}

#endif

