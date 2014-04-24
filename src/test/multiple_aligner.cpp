/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "multiple_aligner.hpp"
#include "PairAligner.hpp"
#include "global.hpp"

using namespace bloomrepeats;

BOOST_AUTO_TEST_CASE (multiple_aligner_eq) {
    Strings seqs((3));
    seqs[0] = "ATTT";
    seqs[1] = "ATTT";
    seqs[2] = "ATTT";
    PairAligner pa(-1, 15, 1);
    multiple_aligner(seqs, &pa);
    BOOST_CHECK(seqs[0] == "ATTT");
    BOOST_CHECK(seqs[1] == "ATTT");
    BOOST_CHECK(seqs[2] == "ATTT");
}

BOOST_AUTO_TEST_CASE (multiple_aligner_n) {
    Strings seqs((3));
    seqs[0] = "ATTT";
    seqs[1] = "ANTT";
    seqs[2] = "ATTT";
    PairAligner pa(-1, 15, 1);
    multiple_aligner(seqs, &pa);
    BOOST_CHECK(seqs[0] == "ATTT");
    BOOST_CHECK(seqs[1] == "ANTT");
    BOOST_CHECK(seqs[2] == "ATTT");
}

BOOST_AUTO_TEST_CASE (multiple_aligner_gap) {
    Strings seqs((3));
    seqs[0] = "ATGC";
    seqs[1] = "AGC";
    seqs[2] = "ATGC";
    PairAligner pa(-1, 15, 1);
    multiple_aligner(seqs, &pa);
    BOOST_CHECK(seqs[0] == "ATGC");
    BOOST_CHECK(seqs[1] == "A-GC");
    BOOST_CHECK(seqs[2] == "ATGC");
}

BOOST_AUTO_TEST_CASE (multiple_aligner_gap_2) {
    Strings seqs((3));
    seqs[0] = "ATGC";
    seqs[1] = "AGC";
    seqs[2] = "ATG";
    PairAligner pa(-1, 15, 1);
    multiple_aligner(seqs, &pa);
    BOOST_CHECK(seqs[0] == "ATGC");
    BOOST_CHECK(seqs[1] == "A-GC");
    BOOST_CHECK(seqs[2] == "ATG-");
}

BOOST_AUTO_TEST_CASE (multiple_aligner_gap_3) {
    Strings seqs((3));
    seqs[0] = "TGC";
    seqs[1] = "CTG";
    seqs[2] = "ATG";
    PairAligner pa(-1, 15, 1);
    multiple_aligner(seqs, &pa);
    BOOST_CHECK(seqs[0] == "-TGC");
    BOOST_CHECK(seqs[1] == "CTG-");
    BOOST_CHECK(seqs[2] == "ATG-");
}

