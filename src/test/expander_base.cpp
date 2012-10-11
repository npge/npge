/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "PairAligner.hpp"
#include "ExpanderBase.hpp"

BOOST_AUTO_TEST_CASE (ExpanderBase_aligned) {
    using namespace bloomrepeats;
    ExpanderBase base;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tggtccga|tggtccga");
    base.aligner().set_max_errors(0); // eq
    base.set_batch(100);
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 0), Fragment(s1, 8, 8)));
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 0), Fragment(s1, 3, 3)));
    BOOST_CHECK(!base.aligned(Fragment(s1, 0, 0), Fragment(s1, 1, 1)));
    base.set_batch(1);
    BOOST_CHECK(!base.aligned(Fragment(s1, 0, 0), Fragment(s1, 1, 1)));
    base.set_batch(100);
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 0), Fragment(s1, 0, 0)));
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 0), Fragment(s1, 7, 7, -1)));
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 7), Fragment(s1, 8, 15)));
    base.set_batch(3);
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 7), Fragment(s1, 8, 15)));
    base.set_batch(1);
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 7), Fragment(s1, 8, 15)));
    base.set_batch(100);
    BOOST_CHECK(!base.aligned(Fragment(s1, 0, 7), Fragment(s1, 8, 14)));
    base.set_aligner(PairAligner(1, 1, 1));
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 7), Fragment(s1, 8, 14)));
    base.set_batch(3);
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 7), Fragment(s1, 8, 14)));
    base.set_batch(1);
    BOOST_CHECK(base.aligned(Fragment(s1, 0, 7), Fragment(s1, 8, 14)));
}

