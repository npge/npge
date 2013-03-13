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

SeqStorage::SeqStorage(SequenceType seq_type):
    seq_type_(seq_type)
{ }

void SeqStorage::add_options_impl(po::options_description& desc) const {
    std::string storage = seq_type() == ASIS_SEQUENCE ? "asis" : "compact";
    add_unique_options(desc)
    ("seq-storage", po::value<std::string>()->default_value(storage),
     "way of storing sequence in memory ('asis' or 'compact')");
   ;
}

void SeqStorage::apply_options_impl(const po::variables_map& vm) {
    std::string storage = vm[prefixed("seq-storage")].as<std::string>();
    if (storage == "asis") {
        set_seq_type(ASIS_SEQUENCE);
    } else if (storage == "compact") {
        set_seq_type(COMPACT_SEQUENCE);
    } else {
        throw Exception("'" + prefixed("seq-storage") + "'"
                        "must be 'asis' or 'compact'");
    }
}

SequencePtr SeqStorage::create_sequence() const {
    return Sequence::new_sequence(seq_type());
}

}

