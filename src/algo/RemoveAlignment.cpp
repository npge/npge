/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "RemoveAlignment.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

namespace bloomrepeats {

void RemoveAlignment::run_impl() const {
    BOOST_FOREACH (Block* b, *block_set()) {
        BOOST_FOREACH (Fragment* f, *b) {
            if (f->row()) {
                f->set_row(0);
            }
        }
    }
}

}

