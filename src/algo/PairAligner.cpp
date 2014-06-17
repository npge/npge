/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <utility>
#include <boost/foreach.hpp>
#include "boost-xtime.hpp"
#include <boost/thread/tss.hpp>

#include "PairAligner.hpp"
#include "GeneralAligner.hpp"
#include "Meta.hpp"
#include "tss_meta.hpp"

namespace npge {

struct PairAlignerContents {
    const char* first_start_;
    const char* second_start_;
    int first_size_, second_size_;
    int mismatch_;

    PairAlignerContents():
        first_start_(0), second_start_(0),
        first_size_(0), second_size_(0),
        mismatch_(1) {
    }

    int first_size() const {
        return first_size_;
    }

    int second_size() const {
        return second_size_;
    }

    int substitution(int row, int col) const {
        return (first_start_[row] == second_start_[col] &&
                first_start_[row] != 'N') ? 0 : mismatch_;
    }
};

struct PairAligner::Impl {
    PairAlignerContents pac_;
    GeneralAligner<PairAlignerContents> ga_;
    bool no_tail_;
};

PairAligner::PairAligner(Meta* meta):
    impl_(new Impl) {
    Meta* m = meta ?: tss_meta();
    GeneralAligner<PairAlignerContents>& g = impl_->ga_;
    g.set_max_errors(m->get_opt("ALIGNER_MAX_ERRORS").as<int>());
    g.set_gap_range(m->get_opt("ALIGNER_GAP_RANGE").as<int>());
    g.set_gap_penalty(m->get_opt("ALIGNER_GAP_PENALTY").as<int>());
    impl_->no_tail_ = true;
}

PairAligner::~PairAligner() {
    delete impl_;
}

boost::thread_specific_ptr<PairAligner> local_aligner_;

PairAligner* PairAligner::default_aligner() {
    if (local_aligner_.get() == 0) {
        local_aligner_.reset(new PairAligner());
    }
    return local_aligner_.get();
}

void PairAligner::set_first(const char* start, int size) {
    impl_->pac_.first_start_ = start;
    impl_->pac_.first_size_ = size;
}

void PairAligner::set_second(const char* start, int size) {
    impl_->pac_.second_start_ = start;
    impl_->pac_.second_size_ = size;
}

int PairAligner::gap_range() const {
    return impl_->ga_.gap_range();
}

void PairAligner::set_gap_range(int gap_range) {
    impl_->ga_.set_gap_range(gap_range);
}

int PairAligner::max_errors() const {
    return impl_->ga_.max_errors();
}

void PairAligner::set_max_errors(int max_errors) {
    impl_->ga_.set_max_errors(max_errors);
}

int PairAligner::gap_penalty() const {
    return impl_->ga_.gap_penalty();
}

void PairAligner::set_gap_penalty(int gap_penalty) {
    return impl_->ga_.set_gap_penalty(gap_penalty);
}

int PairAligner::mismatch_penalty() const {
    return impl_->pac_.mismatch_;
}

void PairAligner::set_mismatch_penalty(int mismatch_penalty) {
    impl_->pac_.mismatch_ = mismatch_penalty;
}

bool PairAligner::no_tail() const {
    return impl_->no_tail_;
}

void PairAligner::set_no_tail(bool no_tail) {
    impl_->no_tail_ = no_tail;
}

bool PairAligner::local() const {
    return impl_->ga_.local();
}

void PairAligner::set_local(bool local) {
    impl_->ga_.set_local(local);
}

void PairAligner::align(int& first_last, int& second_last,
                        std::string* first_str,
                        std::string* second_str,
                        PairAlignment* alignment,
                        char gap) const {
    impl_->ga_.set_contents(impl_->pac_);
    impl_->ga_.align(first_last, second_last);
    if (no_tail()) {
        impl_->ga_.cut_tail(first_last, second_last);
    }
    if (first_str || second_str || alignment) {
        export_alignment(first_last, second_last,
                         first_str, second_str, alignment, gap);
    }
}

bool PairAligner::aligned(const std::string& first,
                          const std::string& second,
                          int* fl, int* sl) {
    set_first(first.c_str(), first.size());
    set_second(second.c_str(), second.size());
    int first_last, second_last;
    bool old_no_tail =  no_tail();
    set_no_tail(false);
    align(first_last, second_last);
    set_no_tail(old_no_tail);
    bool result = first_last == first.size() - 1;
    if (impl_->ga_.in(first.size() - 1, second.size() - 1)) {
        int e = impl_->ga_.at(first.size() - 1, second.size() - 1);
        bool ok = (e <= impl_->ga_.max_errors());
        result &= ok;
    }
    if (fl && sl) {
        if (no_tail()) {
            impl_->ga_.cut_tail(first_last, second_last);
        }
        *fl = first_last;
        *sl = second_last;
    }
    return result;
}

void PairAligner::export_alignment(int row, int col,
                                   std::string* first_str,
                                   std::string* second_str,
                                   PairAlignment* alignment,
                                   char gap) const {
    PairAlignment temp_alignment;
    if (!alignment) {
        alignment = &temp_alignment;
    }
    impl_->ga_.export_alignment(row, col, *alignment);
    if (first_str && second_str) {
        std::string& s1 = *first_str;
        std::string& s2 = *second_str;
        const PairAlignerContents& pc = impl_->pac_;
        typedef std::pair<int, int> Match;
        BOOST_FOREACH (const Match& m, *alignment) {
            s1 += (m.first == -1) ? gap : pc.first_start_[m.first];
            s2 += (m.second == -1) ? gap : pc.second_start_[m.second];
        }
    }
}

}

