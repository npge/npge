/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_CONFIG_HPP_
#define BR_CONFIG_HPP_

namespace bloomrepeats {

const int MIN_LENGTH = 100;
const double MIN_IDENTITY = 0.9;
const double MAX_EVALUE = 0.001;
const double MAX_SPREADING = 0.2;
const double MAX_GAPS = 0.2;
const int EXPANDER_BATCH = 100;
const int EXPANDER_MAX_OVERLAP = 200;
const int ALIGNER_MAX_ERRORS = EXPANDER_BATCH* (1 - MIN_IDENTITY);
const int ALIGNER_GAP_RANGE = 5;
const int ALIGNER_GAP_PENALTY = 2;
const int BOUNDARIES_MIN_DISTANCE = 100;
const int JOINER_MAX_DIST = 100;
const double JOINER_RATIO_TO_FRAGMENT = 0.5;
const double JOINER_GAP_RATIO = 1.5;

}

#endif

