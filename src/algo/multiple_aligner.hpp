/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_MULTIPLE_ALIGNER_HPP_
#define BR_MULTIPLE_ALIGNER_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Align multiple sequences */
void multiple_aligner(Strings& seqs, PairAligner* pa);

}

#endif

