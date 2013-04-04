/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CUT_GAPS_HPP_
#define BR_CUT_GAPS_HPP_

#include "BlocksJobs.hpp"
#include "RowStorage.hpp"

namespace bloomrepeats {

/** Cut longest terminal gap.

Alignment is preserved.

Example:
Before: "-aaaaa-----a-". After: "aaaaa-----a".

If a fragment consists only of gaps after cut,
it is removed from block.
*/
class CutGaps : public BlocksJobs, public RowStorage {
public:
    /** The mode */
    enum Mode {
        /** Cut all columns from the beginning until first gapless column.
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
        */
        STRICT,

        /** Find end gap of max length and cut those columns.
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
        PERMISSIVE
    };

    /** Constructor */
    CutGaps(Mode mode = PERMISSIVE);

    /** Get mode */
    Mode mode() const {
        return mode_;
    };

    /** Set mode */
    void set_mode(Mode mode) {
        mode_ = mode;
    };

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool apply_to_block_impl(Block* block) const;

    const char* name_impl() const;

private:
    Mode mode_;
};

}

#endif

