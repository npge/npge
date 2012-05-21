/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <boost/test/unit_test.hpp>

#include "make_hash.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"

BOOST_AUTO_TEST_CASE (hash_main) {
    using namespace bloomrepeats;
    std::string s("CGCAtacccTGCGgcaGGGTcaGGGC");
    Sequence::to_atgc(s);
    BOOST_REQUIRE(make_hash(s.c_str(), 4) == make_hash(s.c_str() + 12, 4, -1));
    BOOST_WARN(make_hash(s.c_str(), 4) != make_hash(s.c_str() + 12, 4, 1));
    BOOST_WARN(make_hash(s.c_str() + 16, 4) != make_hash(s.c_str() + 22, 4));
}

BOOST_AUTO_TEST_CASE (hash_fragment) {
    using namespace bloomrepeats;
    std::string s("cgcataccctgcggcagggtcagggc");
    Sequence::to_atgc(s);
    SequencePtr s1 = boost::make_shared<InMemorySequence>(s);
    BOOST_REQUIRE(Fragment(s1, 0, 3).hash() == make_hash(s.c_str(), 4));
    BOOST_REQUIRE(Fragment(s1, 1, 4).hash() == make_hash(s.c_str() + 1, 4));
    BOOST_REQUIRE(Fragment(s1, 0, 3, -1).hash() ==
                  make_hash(s.c_str() + 3, 4, -1));
}

