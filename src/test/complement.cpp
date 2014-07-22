/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>

#include "complement.hpp"
#include "make_hash.hpp"
#include "global.hpp"

BOOST_AUTO_TEST_CASE (complement_char) {
    using namespace npge;
    BOOST_CHECK(complement('A') == 'T');
    BOOST_CHECK(complement('T') == 'A');
    BOOST_CHECK(complement('G') == 'C');
    BOOST_CHECK(complement('C') == 'G');
}

BOOST_AUTO_TEST_CASE (complement_string) {
    using namespace npge;
    std::string data;
    data = "AAA";
    complement(data);
    BOOST_CHECK(data == "TTT");
    //
    data = "ATGC";
    complement(data);
    BOOST_CHECK(data == "GCAT");
    //
    data = "ACGT";
    complement(data);
    BOOST_CHECK(data == "ACGT");
    //
    data = "A-CG-T";
    complement(data);
    BOOST_CHECK(data == "A-CG-T");
    //
    data = "TT~~A";
    complement(data);
    BOOST_CHECK(data == "T~~AA");
}

static npge::hash_t make_hash(const std::string& s) {
    using namespace npge;
    return make_hash(s.c_str(), s.length());
}

BOOST_AUTO_TEST_CASE (complement_hash_test) {
    using namespace npge;
    BOOST_CHECK(make_hash("AAA") == make_hash("AAA"));
    BOOST_CHECK(complement_hash(make_hash("AAA"), 3) ==
                make_hash("TTT"));
    BOOST_CHECK(complement_hash(make_hash("ATGC"), 4) ==
                make_hash("GCAT"));
    BOOST_CHECK(complement_hash(make_hash("ATGCAAAAG"), 9) ==
                make_hash("CTTTTGCAT"));
    Strings tests = boost::assign::list_of
                    ("CTTTTGCAT")
                    ("CATATCGGTAAGGTGACATTA")
                    ("CGCCCGCGAGTGGAATACAACACCTCCACCTG")
                    ("CGCCCGCGAGTGGAATACAACACCTCCACCT")
                    ("CGCCCGCGAGTGGAATACAACACCTCCACC")
                    ("CGCCCGCGAGTGGAATACAACACCTC");
    BOOST_FOREACH (std::string t, tests) {
        std::string t_cmpl = t;
        complement(t_cmpl);
        BOOST_CHECK(complement_hash(make_hash(t), t.size()) ==
                    make_hash(t_cmpl));
    }
}

