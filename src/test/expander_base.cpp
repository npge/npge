/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "PairAligner.hpp"
#include "ExpanderBase.hpp"
#include "Pipe.hpp"

BOOST_AUTO_TEST_CASE (ExpanderBase_aligned) {
    using namespace npge;
    Pipe _;
    add_expander_options(&_);
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccga|tggtccga");
    _.set_opt_value("max-errors", 0); // eq
    _.set_opt_value("batch", 100);
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 0), Fragment(s1, 8, 8)));
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 0), Fragment(s1, 3, 3)));
    BOOST_CHECK(!aligned(&_, Fragment(s1, 0, 0), Fragment(s1, 1, 1)));
    _.set_opt_value("batch", 1);
    BOOST_CHECK(!aligned(&_, Fragment(s1, 0, 0), Fragment(s1, 1, 1)));
    _.set_opt_value("batch", 100);
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 0), Fragment(s1, 0, 0)));
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 0), Fragment(s1, 7, 7, -1)));
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 7), Fragment(s1, 8, 15)));
    _.set_opt_value("batch", 3);
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 7), Fragment(s1, 8, 15)));
    _.set_opt_value("batch", 1);
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 7), Fragment(s1, 8, 15)));
    _.set_opt_value("batch", 100);
    BOOST_CHECK(!aligned(&_, Fragment(s1, 0, 7), Fragment(s1, 8, 14)));
    _.set_opt_value("max-errors", 1);
    _.set_opt_value("gap-penalty", 1);
    _.set_opt_value("gap-range", 1);
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 7), Fragment(s1, 8, 14)));
    _.set_opt_value("batch", 3);
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 7), Fragment(s1, 8, 14)));
    _.set_opt_value("batch", 1);
    BOOST_CHECK(aligned(&_, Fragment(s1, 0, 7), Fragment(s1, 8, 14)));
}

