/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_SEQ_STORAGE_HPP_
#define BR_SEQ_STORAGE_HPP_

#include "OptionsPrefix.hpp"

namespace bloomrepeats {

/** Add sequence storage configuration to a processor */
void add_seq_storage_options(Processor* processor);

/** Return sequence type processor uses */
SequenceType seq_type(const Processor* processor);

/** Create sequence using sequence storage configuration options */
SequencePtr create_sequence(const Processor* processor);

}

#endif

