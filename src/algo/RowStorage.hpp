/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ROW_STORAGE_HPP_
#define BR_ROW_STORAGE_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Add row storage configuration to a processor */
void add_row_storage_options(Processor* processor);

/** Return row type processor uses */
RowType row_type(const Processor* processor);

/** Create row using row storage configuration options */
AlignmentRow* create_row(const Processor* processor);

}

#endif

