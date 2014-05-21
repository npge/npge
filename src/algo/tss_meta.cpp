/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include "boost-xtime.hpp"
#include <boost/thread/tss.hpp>

#include "tss_meta.hpp"
#include "Meta.hpp"

namespace bloomrepeats {

static boost::thread_specific_ptr<Meta> tss_meta_;

Meta* tss_meta() {
    if (tss_meta_.get() == 0) {
        void* ptr = malloc(sizeof(Meta));
        Meta* meta = reinterpret_cast<Meta*>(ptr);
        tss_meta_.reset(meta);
        new(meta) Meta;
    }
    return tss_meta_.get();
}

AnyAs tss_go(const std::string& key, const AnyAs& dflt) {
    return tss_meta()->get_opt(key, dflt);
}

}

