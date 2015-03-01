/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Decimal.hpp"

BOOST_AUTO_TEST_CASE (Decimal_main) {
    using namespace npge;
    BOOST_CHECK(D(1.01) + Decimal() == D(1.01));
    BOOST_CHECK(D(1.01) + Decimal() >= D(1.01));
    BOOST_CHECK(D(1.01) + Decimal() <= D(1.01));
    BOOST_CHECK(D(1.01) > D(1.001));
    BOOST_CHECK(D(1.001) < D(1.01));
    BOOST_CHECK(D(-1.001) < D(1.01));
    BOOST_CHECK(D(-1.001) > D(-1.01));
    BOOST_CHECK(D(0.001) > Decimal());
    BOOST_CHECK(D(0.001).to_s() == "0.001" ||
                D(0.001).to_s() == ".001");
    BOOST_CHECK(D(0).to_s() == "0.0" ||
                D(0).to_s() == "0." ||
                D(0).to_s() == ".0" ||
                D(0).to_s() == "0");
    BOOST_CHECK(Decimal(Decimal().to_s()) == Decimal());
    BOOST_CHECK(D(0.0) == Decimal());
    BOOST_CHECK(D(.0) == Decimal());
    BOOST_CHECK(D(0.) == Decimal());
    BOOST_CHECK(D(0) == Decimal());
    BOOST_CHECK(D(1) / 2 == D(0.5));
    BOOST_CHECK(D(0.75) * D(0.5) == D(0.375));
    int tenth = Decimal::sub_point / 10;
    BOOST_CHECK(D(1.3).fraction() == 3 * tenth);
    BOOST_CHECK(D(-1.3).fraction() == 3 * tenth);
    BOOST_CHECK(D(1.3).to_i() == 1);
    BOOST_CHECK(D(-1.3).to_i() == -1);
    BOOST_CHECK(D(1.3).round() == 1);
    BOOST_CHECK(D(-1.3).round() == -1);
    BOOST_CHECK(D(1.7).fraction() == 7 * tenth);
    BOOST_CHECK(D(-1.7).fraction() == 7 * tenth);
    BOOST_CHECK(D(1.7).to_i() == 1);
    BOOST_CHECK(D(-1.7).to_i() == -1);
    BOOST_CHECK(D(1.7).round() == 2);
    BOOST_CHECK(D(-1.7).round() == -2);
    BOOST_CHECK(D(0.9) * 100 == 90);
    BOOST_CHECK((D(0.9) * 100).to_i() == 90);
}

