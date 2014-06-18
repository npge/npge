/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "SimilarAligner.hpp"
#include "DummyAligner.hpp"

BOOST_AUTO_TEST_CASE (Aligner_test) {
    using namespace npge;
    SimilarAligner sa;
    BOOST_CHECK(sa.test(/* gaps */ false));
    BOOST_CHECK(sa.test(/* gaps */ true));
    DummyAligner da;
    BOOST_CHECK(da.test());
    BOOST_CHECK(da.test(/* gaps */ false));
    BOOST_CHECK(!da.test(/* gaps */ true));
}

