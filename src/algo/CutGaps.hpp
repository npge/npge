/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CUT_GAPS_HPP_
#define BR_CUT_GAPS_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Cut longest terminal gap.

Alignment is preserved.

Example:
Before: "-aaaaa-----a-". After: "aaaaa-----a".

If a fragment consists only of gaps after cut,
it is removed from block.

<h4>Modes</h4>
<h5>strict</h5>
Cut all columns from the beginning until first gapless column.
If there is no gapless column, then the block is cleared.
Example:
\code
Before:
    ----aaaa----
    --aaaaaaa---
    --aa-aaaa---

After:
         aaa
         aaa
         aaa
\endcode

<h5>permissive</h5>
Find end gap of max length and cut those columns.
If end gap of max length from begin and from end overlap,
block is cleared.

Example:
\code
Before:
    ----aaaa----
    --aaaaaaa---
    --aa-aaaa---

After:
        aaaa
        aaaa
        -aaa
\endcode
*/
class CutGaps : public BlocksJobs {
public:
    /** Constructor */
    CutGaps(bool strict = false);

    /** Do the job and return if the block was changed */
    bool cut_gaps(Block* block) const;

protected:
    void process_block_impl(Block* block, ThreadData*) const;

    const char* name_impl() const;
};

}

#endif

