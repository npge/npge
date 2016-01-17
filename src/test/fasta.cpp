/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/test/unit_test.hpp>

#include "BlockSet.hpp"
#include "Sequence.hpp"

BOOST_AUTO_TEST_CASE (fasta_main) {
    using namespace npge;
    std::stringstream ss(">name description\nAAA");
    BlockSet bs;
    ss >> bs;
    BOOST_REQUIRE(bs.seqs().size() == 1);
    BOOST_CHECK(bs.seqs()[0]->name() == "name");
    BOOST_CHECK(bs.seqs()[0]->description() == "description");
}

BOOST_AUTO_TEST_CASE (fasta_tab) {
    using namespace npge;
    std::stringstream ss(">name\tdescription\nAAA");
    BlockSet bs;
    ss >> bs;
    BOOST_REQUIRE(bs.seqs().size() == 1);
    BOOST_CHECK(bs.seqs()[0]->name() == "name");
    BOOST_CHECK(bs.seqs()[0]->description() == "description");
}

BOOST_AUTO_TEST_CASE (fasta_2spaces) {
    using namespace npge;
    std::stringstream ss(">name "
                         " description\nAAA");
    BlockSet bs;
    ss >> bs;
    BOOST_REQUIRE(bs.seqs().size() == 1);
    BOOST_CHECK(bs.seqs()[0]->name() == "name");
    BOOST_CHECK(bs.seqs()[0]->description() == "description");
}

BOOST_AUTO_TEST_CASE (fasta_words) {
    using namespace npge;
    std::stringstream ss(">name \t a b\tc\nAAA");
    BlockSet bs;
    ss >> bs;
    BOOST_REQUIRE(bs.seqs().size() == 1);
    BOOST_CHECK(bs.seqs()[0]->name() == "name");
    BOOST_CHECK(bs.seqs()[0]->description() == "a b\tc");
}

