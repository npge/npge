/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/bind.hpp>

#include "MetaAligner.hpp"
#include "ExternalAligner.hpp"
#include "MultipleAligner.hpp"
#include "SimilarAligner.hpp"
#include "DummyAligner.hpp"
#include "throw_assert.hpp"
#include "global.hpp"
#include "config.hpp"

namespace bloomrepeats {

bool MetaAligner::check_type(std::string& m) const {
    std::string a_type = opt_value("aligner-type").as<std::string>();
    if (a_type == "external") {
        aligner_ = external_;
    } else if (a_type == "multiple") {
        aligner_ = multiple_;
    } else if (a_type == "similar") {
        aligner_ = similar_;
    } else if (a_type == "dummy") {
        aligner_ = dummy_;
    } else {
        m = "bad aligner-type: " + a_type;
        return false;
    }
    return true;
}

MetaAligner::MetaAligner() {
    external_ = new ExternalAligner;
    external_->set_parent(this);
    multiple_ = new MultipleAligner;
    multiple_->set_parent(this);
    similar_ = new SimilarAligner;
    similar_->set_parent(this);
    dummy_ = new DummyAligner;
    dummy_->set_parent(this);
    aligner_ = 0;
    add_opt("aligner-type", "Type of aligner "
            "(external, multiple, similar, dummy)",
            std::string("similar"));
    add_opt_check(boost::bind(&MetaAligner::check_type, this, _1));
}

void MetaAligner::align_seqs_impl(Strings& seqs) const {
    if (!aligner_) {
        std::string m;
        ASSERT_TRUE(check_type(m));
    }
    ASSERT_TRUE(aligner_);
    aligner_->align_seqs(seqs);
}

const char* MetaAligner::name_impl() const {
    return "Align blocks with one of possible aligners";
}

}

