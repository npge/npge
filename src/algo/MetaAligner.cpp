/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "MetaAligner.hpp"
#include "ExternalAligner.hpp"
#include "SimilarAligner.hpp"
#include "DummyAligner.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

bool MetaAligner::check_type(std::string& m) const {
    std::string a_type = opt_value("aligner-type").as<std::string>();
    if (a_type == last_aligners_) {
        return true;
    }
    using namespace boost::algorithm;
    Strings aligners;
    split(aligners, a_type, is_any_of(","));
    ASSERT_GTE(aligners.size(), 1);
    BOOST_FOREACH (const std::string& aligner, aligners) {
        aligner_ = 0;
        BOOST_FOREACH (AbstractAligner* a, aligners_) {
            if (aligner == a->aligner_type()) {
                aligner_ = a;
                break;
            }
        }
        if (!aligner_) {
            m = "bad aligner-type: " + aligner;
            return false;
        }
        if (aligners.size() == 1 || aligner_->test()) {
            break;
        }
    }
    ASSERT_TRUE(aligner_);
    if (aligners.size() >= 2) {
        write_log("Selected aligner: " +
                  aligner_->aligner_type());
    }
    last_aligners_ = a_type;
    return true;
}

MetaAligner::MetaAligner() {
    add_aligner(new MafftAligner);
    add_aligner(new MuscleAligner);
    add_aligner(new SimilarAligner);
    add_aligner(new DummyAligner);
    aligner_ = 0;
    add_gopt("aligner-type", "Type of aligner "
             "(external, mafft, muscle, "
             "similar, dummy). Specify several types, "
             "separated by comma, the first working one "
             "will be used or the last one if all fail.",
             "ALIGNER");
    add_opt_check(boost::bind(&MetaAligner::check_type, this, _1));
}

void MetaAligner::add_aligner(AbstractAligner* aligner) {
    aligner->set_parent(this);
    aligners_.push_back(aligner);
}

void MetaAligner::align_seqs_impl(Strings& seqs) const {
    if (!aligner_) {
        std::string m;
        bool ok = check_type(m);
        ASSERT_TRUE(ok);
    }
    ASSERT_TRUE(aligner_);
    aligner_->align_seqs(seqs);
}

std::string MetaAligner::aligner_type() const {
    return "meta";
}

const char* MetaAligner::name_impl() const {
    return "Align blocks with one of possible aligners";
}

}

