/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "InternalAligner.hpp"
#include "ExpanderBase.hpp"
#include "PairAligner.hpp"
#include "multiple_aligner.hpp"

namespace bloomrepeats {

InternalAligner::InternalAligner() {
    add_expander_options(this);
}

void InternalAligner::align_seqs_impl(Strings& seqs) const {
    PairAligner pa;
    apply_pair_aligner_options(&pa, this);
    multiple_aligner(seqs, &pa);
}

const char* InternalAligner::name_impl() const {
    return "Internal aligner";
}

}

