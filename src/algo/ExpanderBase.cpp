/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <climits>

#include "ExpanderBase.hpp"
#include "PairAligner.hpp"
#include "Fragment.hpp"
#include "Processor.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"
#include "config.hpp"

namespace bloomrepeats {

void add_expander_options(Processor* p) {
    p->add_opt("batch", "batch size for pair aligner", EXPANDER_BATCH);
    p->add_opt("gap-range", "Max distance from main diagonal of "
               "considered states of pair alignment. "
               "The more gap_range, the more time.", ALIGNER_GAP_RANGE);
    p->add_opt("max-errors", "Max number of errors in pair alignment",
               ALIGNER_MAX_ERRORS);
    p->add_opt("gap-penalty", "Gap open or extension penalty",
               ALIGNER_GAP_PENALTY);
    p->add_opt("mismatch-penalty", "Mismatch penalty",
               ALIGNER_MISMATCH_PENALTY);
}

void apply_pair_aligner_options(PairAligner* pa, const Processor* p) {
    pa->set_max_errors(p->opt_value("max-errors").as<int>());
    pa->set_gap_range(p->opt_value("gap-range").as<int>());
    pa->set_gap_penalty(p->opt_value("gap-penalty").as<int>());
}

bool aligned(const Processor* p,
             const Fragment& f1, const Fragment& f2) {
    int batch = p->opt_value("batch").as<int>();
    int max_errors = p->opt_value("max-errors").as<int>();
    int gap_range = p->opt_value("gap-range").as<int>();
    int gap_penalty = p->opt_value("gap-penalty").as<int>();
    PairAligner pa(max_errors, gap_range, gap_penalty);
    int f1_last = -1, f2_last = -1;
    ASSERT_LTE(f1.length(), INT_MAX);
    ASSERT_LTE(f2.length(), INT_MAX);
    while (f1_last < int(f1.length()) - 1 &&
            f2_last < int(f2.length()) - 1) {
        int sub_this_last, sub_other_last;
        int this_min = std::min(int(f1.length()) - 1, f1_last + 1);
        int this_max = std::min(int(f1.length()) - 1, f1_last + batch);
        int other_min = std::min(int(f2.length()) - 1, f2_last + 1);
        int other_max = std::min(int(f2.length()) - 1, f2_last + batch);
        if (!pa.aligned(f1.substr(this_min, this_max),
                        f2.substr(other_min, other_max),
                        &sub_this_last, &sub_other_last)) {
            return false;
        }
        if (sub_this_last == -1 || sub_other_last == -1) {
            return false;
        }
        f1_last += sub_this_last + 1;
        f2_last += sub_other_last + 1;
    }
    ASSERT_LT(f1_last, int(f1.length()));
    ASSERT_LT(f2_last, int(f2.length()));
    return pa.aligned(f1.substr(f1_last, f1.length() - 1),
                      f2.substr(f2_last, f2.length() - 1));
}

}

