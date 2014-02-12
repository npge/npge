/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/bind.hpp>

#include "SeqStorage.hpp"
#include "Processor.hpp"
#include "Sequence.hpp"

namespace bloomrepeats {

static bool check_seq_type(std::string& message, Processor* p) {
    std::string seq_type = p->opt_value("seq-storage").as<std::string>();
    if (seq_type != "asis" && seq_type != "compact") {
        message = "seq-storage must be 'asis' or 'compact'";
        return false;
    }
    return true;
}

void add_seq_storage_options(Processor* p) {
    p->add_opt("seq-storage",
               "way of storing sequences in memory ('asis' or 'compact')",
               std::string("compact"));
    p->add_opt_check(boost::bind(check_seq_type, _1, p));
}

SequenceType seq_type(const Processor* p) {
    return (p->opt_value("seq-storage").as<std::string>() == "asis") ?
        ASIS_SEQUENCE : COMPACT_SEQUENCE;
}

SequencePtr create_sequence(const Processor* p) {
    return Sequence::new_sequence(seq_type(p));
}

}

