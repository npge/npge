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
    Sequence::to_atgcn(s);
    BOOST_CHECK(make_hash(s.c_str(), 4) == make_hash(s.c_str() + 12, 4, -1));
    BOOST_WARN(make_hash(s.c_str(), 4) != make_hash(s.c_str() + 12, 4, 1));
    BOOST_WARN(make_hash(s.c_str() + 16, 4) != make_hash(s.c_str() + 22, 4));
}

BOOST_AUTO_TEST_CASE (hash_fragment) {
    using namespace bloomrepeats;
    std::string s("CGCATACCCTGCGGCAGGGTCAGGGC");
    Sequence::to_atgcn(s);
    SequencePtr s1 = boost::make_shared<InMemorySequence>(s);
    BOOST_CHECK(Fragment(s1, 0, 3).hash() == make_hash(s.c_str(), 4));
    BOOST_CHECK(Fragment(s1, 1, 4).hash() == make_hash(s.c_str() + 1, 4));
    BOOST_CHECK(Fragment(s1, 0, 3, -1).hash() ==
                make_hash(s.c_str() + 3, 4, -1));
}

BOOST_AUTO_TEST_CASE (hash_reuse_hash) {
    using namespace bloomrepeats;
    std::string s("CGCATACCCTGCGGCAGGGTCAGGGC");
    Sequence::to_atgcn(s);
    size_t h = make_hash(s.c_str(), 4);
    BOOST_CHECK(reuse_hash(h, 4, s[0], s[4]) == make_hash(s.c_str() + 1, 4));
}

BOOST_AUTO_TEST_CASE (hash_reuse_hash_full) {
    using namespace bloomrepeats;
    std::string s("GATCCTCGATTAACAGTTTGGCCTGTTCCTATGTATGCCCTACTCCAAATGGT"
                  "GCCAACTGGATCAATCCTCAGTGCCGCGGGAATCATGTCTTTATTTATGCTTT"
                  "TCAGCTCTGCGAACTTAGGCTCAGCACAAGATTTAAGCGAGAAGCGAAAGCTG"
                  "ACCGGCAGGGGGGGCACGGTTAATAACTAAGACTGTAGCGTGACAAACGGACC");
    SequencePtr s1 = boost::make_shared<InMemorySequence>(s);
    for (int i = 70; i < 100; i++) {
        for (int length = 1; length < 100; length++) {
            for (int fr_ori = -1; fr_ori <= 1; fr_ori += 2) {
                for (int move_ori = -1; move_ori <= 1; move_ori += 2) {
                    Fragment f(s1, i, i + length - 1, fr_ori);
                    size_t h = f.hash();
                    bool forward = move_ori == fr_ori;
                    char remove = f.at(forward ? 0 : -1);
                    char add = f.raw_at(forward ? length : -1);
                    size_t reused = reuse_hash(h, length, remove, add, forward);
                    f.set_min_pos(f.min_pos() + move_ori);
                    f.set_max_pos(f.max_pos() + move_ori);
                    size_t new_h = f.hash();
                    BOOST_CHECK(reused == new_h);
                }
            }
        }
    }
}

