/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "boost-xtime.hpp"
#include <boost/thread/tss.hpp>

#include "tss_meta.hpp"
#include "Meta.hpp"

namespace npge {

static void do_nothing(Meta*) {
}

static boost::thread_specific_ptr<Meta> tss_meta_(do_nothing);

TssMetaHolder::TssMetaHolder(Meta* meta) {
    prev_ = tss_meta();
    tss_meta_.reset(meta);
}

TssMetaHolder::~TssMetaHolder() {
    tss_meta_.reset(prev_);
}

Meta* tss_meta() {
    return tss_meta_.get();
}

}

