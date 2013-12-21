/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_RAND_NAME_HPP_
#define BR_RAND_NAME_HPP_

#include <string>

namespace bloomrepeats {

/** Return random name */
std::string rand_name(int size);

}

#endif

