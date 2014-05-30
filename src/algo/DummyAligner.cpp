/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <algorithm>
#include <boost/foreach.hpp>

#include "DummyAligner.hpp"

namespace npge {

DummyAligner::DummyAligner() {
}

void DummyAligner::align_seqs_impl(Strings& seqs) const {
    int max_length = seqs.front().length();
    BOOST_FOREACH (const std::string& seq, seqs) {
        max_length = std::max(max_length, int(seq.length()));
    }
    BOOST_FOREACH (std::string& seq, seqs) {
        seq.resize(max_length, '-');
    }
}

const char* DummyAligner::name_impl() const {
    return "Align by adding gaps to sequences";
}

}

