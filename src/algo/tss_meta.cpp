/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/thread/tss.hpp>

#include "tss_meta.hpp"
#include "Meta.hpp"

namespace bloomrepeats {

static boost::thread_specific_ptr<Meta> tss_meta_;

Meta* tss_meta() {
    if (tss_meta_.get() == 0) {
        tss_meta_.reset(new Meta);
    }
    return tss_meta_.get();
}

}

