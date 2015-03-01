/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "SimilarAligner.hpp"
#include "DummyAligner.hpp"
#include "ExternalAligner.hpp"

BOOST_AUTO_TEST_CASE (Aligner_test) {
    using namespace npge;
    SimilarAligner sa;
    BOOST_CHECK(sa.test(/* gaps */ false));
    BOOST_CHECK(sa.test(/* gaps */ true));
    DummyAligner da;
    BOOST_CHECK(da.test());
    BOOST_CHECK(da.test(/* gaps */ false));
    BOOST_CHECK(!da.test(/* gaps */ true));
    ExternalAligner ea_bad;
    // dir command is not an aligner
    ea_bad.set_opt_value("aligner-cmd", std::string("dir"));
    BOOST_CHECK(!ea_bad.test(/* gaps */ false));
    BOOST_CHECK(!ea_bad.test(/* gaps */ true));
}

