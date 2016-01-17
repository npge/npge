/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_UNIQUE_NAMES_HPP_
#define NPGE_UNIQUE_NAMES_HPP_

#include "BlocksJobs.hpp"

namespace npge {

/** Set unique names to all blocks of this blockset.
If name is not default and not unique:
 - "n<num>" is appended with num minimum number to make name unique.

If (name is default or "") and not unique:
 - block_name() is used, if name is null.
 - Block::set_random_name() is called untill the name is unique.

If name of sequece is empty or not unique, it is changed to random.
*/
class UniqueNames : public BlocksJobs {
public:
    /** Constructor */
    UniqueNames();

protected:
    void initialize_work_impl() const;
    void process_block_impl(Block* block, ThreadData* data) const;
    void finish_work_impl() const;
    const char* name_impl() const;

private:
    mutable int genomes_;
};

}

#endif

