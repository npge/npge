/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OVERLAPS_RESOLVER_HPP_
#define BR_OVERLAPS_RESOLVER_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Resolve overlaping fragments.
If some blocks from the block set have overlaping fragments
these two blocks are replaced with one higher (and narrower)
block and several remainder blocks.

Anyway, applying this method guarantees that no blocks of the block set
have overlaping fragments.

Since OverlapsResolver can split blocks,
applying \ref join "join(0)" is recommended.

\verbatim
Input:
    Block 1:
        seq1: ---xxxx----
        seq2: ---xxxx----
    Block 2
        seq2: -----xxxx--
        seq3: -----xxxx--
        seq4: -----xxxx--

Output of OverlapsResolver:
    Block 1:
        seq1: ---xx------
        seq2: ---xx------
    Block 2
        seq2: -------xx--
        seq3: -------xx--
        seq4: -------xx--
    Block 3:
        seq1: -----xx----
        seq2: -----xx----
        seq3: -----xx----
        seq4: -----xx----
\endverbatim
\warning Fragments must be \ref Connector "connected"
   for this to work correctly.
*/
class OverlapsResolver : public Processor {
public:
    /** Return if there are blocks which have overlaping fragments.
    \warning Fragments must be \ref Connector "connected"
    */
    bool overlaps() const;

protected:
    /** Apply the action */
    void run_impl() const;
};

}

#endif

