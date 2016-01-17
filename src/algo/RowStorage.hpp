/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ROW_STORAGE_HPP_
#define NPGE_ROW_STORAGE_HPP_

#include "global.hpp"

namespace npge {

/** Add row storage configuration to a processor */
void add_row_storage_options(Processor* processor);

/** Return row type processor uses */
RowType row_type(const Processor* processor);

/** Create row using row storage configuration options */
AlignmentRow* create_row(const Processor* processor);

}

#endif

