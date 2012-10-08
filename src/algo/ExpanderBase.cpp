/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "ExpanderBase.hpp"

namespace bloomrepeats {

ExpanderBase::ExpanderBase(int batch):
    batch_(batch)
{ }

}

