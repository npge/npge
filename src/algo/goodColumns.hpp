/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_GOOD_COLUMNS_HPP_
#define NPGE_GOOD_COLUMNS_HPP_

#include <vector>

namespace npge {

const int MAX_COLUMN_SCORE = 100;

typedef std::vector<int> Scores;

Scores goodColumns(const char** rows, int nrows, int length,
                   int min_identity, int min_length);

}

#endif
