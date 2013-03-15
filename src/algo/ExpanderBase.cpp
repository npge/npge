/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <climits>

#include "ExpanderBase.hpp"
#include "Fragment.hpp"
#include "po.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

ExpanderBase::ExpanderBase(int batch):
    batch_(batch)
{ }

bool ExpanderBase::aligned(const Fragment& f1, const Fragment& f2) const {
    int f1_last = -1, f2_last = -1;
    BOOST_ASSERT(f1.length() <= INT_MAX);
    BOOST_ASSERT(f2.length() <= INT_MAX);
    while (f1_last < int(f1.length()) - 1 &&
            f2_last < int(f2.length()) - 1) {
        int sub_this_last, sub_other_last;
        int this_min = std::min(int(f1.length()) - 1, f1_last + 1);
        int this_max = std::min(int(f1.length()) - 1, f1_last + batch());
        int other_min = std::min(int(f2.length()) - 1, f2_last + 1);
        int other_max = std::min(int(f2.length()) - 1, f2_last + batch());
        if (!aligner().aligned(f1.substr(this_min, this_max),
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
    BOOST_ASSERT(f1_last < int(f1.length()) && f2_last < int(f2.length()));
    return aligner().aligned(f1.substr(f1_last, f1.length() - 1),
                             f2.substr(f2_last, f2.length() - 1));
}

void ExpanderBase::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("batch", po::value<int>()->default_value(batch()),
     "batch size for pair aligner")
    ("gap-range", po::value<int>()->default_value(aligner().gap_range()),
     "Max distance from main diagonal of considered states of pair alignment. "
     "The more gap_range, the more time.")
    ("max-errors", po::value<int>()->default_value(aligner().max_errors()),
     "Max number of errors in pair alignment")
    ("gap-penalty", po::value<int>()->default_value(aligner().gap_penalty()),
     "Gap open or extension penalty")
   ;
}

void ExpanderBase::apply_options_impl(const po::variables_map& vm) {
    if (vm["batch"].as<int>() < 10) {
        throw Exception("'batch' must be >= 10");
    }
    set_batch(vm["batch"].as<int>());
    aligner().set_gap_range(vm["gap-range"].as<int>());
    aligner().set_max_errors(vm["max-errors"].as<int>());
    aligner().set_gap_penalty(vm["gap-penalty"].as<int>());
}

}

