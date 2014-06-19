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
#include "MultipleAligner.hpp"
#include "SimilarAligner.hpp"
#include "DummyAligner.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

bool MetaAligner::check_type(std::string& m) const {
    std::string a_type = opt_value("aligner-type").as<std::string>();
    using namespace boost::algorithm;
    Strings aligners;
    split(aligners, a_type, is_any_of(","));
    ASSERT_GTE(aligners.size(), 1);
    BOOST_FOREACH (const std::string& aligner, aligners) {
        if (aligner == "external") {
            aligner_ = external_;
        } else if (aligner == "mafft") {
            aligner_ = mafft_;
        } else if (aligner == "muscle") {
            aligner_ = muscle_;
        } else if (aligner == "multiple") {
            aligner_ = multiple_;
        } else if (aligner == "similar") {
            aligner_ = similar_;
        } else if (aligner == "dummy") {
            aligner_ = dummy_;
        } else {
            m = "bad aligner-type: " + aligner;
            return false;
        }
        if (aligner_->test()) {
            break;
        }
    }
    return true;
}

MetaAligner::MetaAligner() {
    external_ = new ExternalAligner;
    external_->set_parent(this);
    mafft_ = new MafftAligner;
    mafft_->set_parent(this);
    muscle_ = new MuscleAligner;
    muscle_->set_parent(this);
    multiple_ = new MultipleAligner;
    multiple_->set_parent(this);
    similar_ = new SimilarAligner;
    similar_->set_parent(this);
    dummy_ = new DummyAligner;
    dummy_->set_parent(this);
    aligner_ = 0;
    add_gopt("aligner-type", "Type of aligner "
             "(external, mafft, muscle, multiple, "
             "similar, dummy). Specify several types, "
             "separated by comma, the first working one "
             "will be used or the last one if all fail.",
             "META_ALIGNER");
    add_opt_check(boost::bind(&MetaAligner::check_type, this, _1));
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

const char* MetaAligner::name_impl() const {
    return "Align blocks with one of possible aligners";
}

}

