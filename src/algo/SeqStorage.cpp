/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "SeqStorage.hpp"
#include "Sequence.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

SeqStorage::SeqStorage(const std::string& storage):
    storage_(storage)
{ }

void SeqStorage::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("seq-storage", po::value<std::string>()->default_value(storage()),
     "way of storing sequence in memory ('asis' or 'compact')");
   ;
}

void SeqStorage::apply_options_impl(const po::variables_map& vm) {
    std::string storage = vm[prefixed("seq-storage")].as<std::string>();
    if (storage != "asis" && storage != "compact") {
        throw Exception("'" + prefixed("seq-storage") + "'"
                        "must be 'asis' or 'compact'");
    }
    set_storage(storage);
}

SequencePtr SeqStorage::create_sequence() const {
    if (storage() == "asis") {
        return SequencePtr(new InMemorySequence);
    } else {
        return SequencePtr(new CompactSequence);
    }
}

}

