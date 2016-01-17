/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/bind.hpp>

#include "SeqStorage.hpp"
#include "Processor.hpp"
#include "Sequence.hpp"

namespace npge {

static bool check_seq_type(std::string& message, Processor* p) {
    std::string st;
    st = p->opt_value("seq-storage").as<std::string>();
    if (st != "asis" && st != "compact" &&
            st != "compact_low_n") {
        message = "seq-storage must be 'asis', 'compact' "
                  "or 'compact_low_n'";
        return false;
    }
    return true;
}

void add_seq_storage_options(Processor* p) {
    p->add_opt("seq-storage",
               "way of storing sequences in memory "
               "('asis', 'compact' or 'compact_low_n')",
               std::string("compact_low_n"));
    p->add_opt_check(boost::bind(check_seq_type, _1, p));
}

SequenceType seq_type(const Processor* p) {
    std::string st;
    st = p->opt_value("seq-storage").as<std::string>();
    return (st == "asis") ? ASIS_SEQUENCE :
           (st == "compact") ? COMPACT_SEQUENCE :
           COMPACT_LOW_N_SEQUENCE;
}

SequencePtr create_sequence(const Processor* p) {
    return Sequence::new_sequence(seq_type(p));
}

}

