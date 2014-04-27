/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "SimilarAligner.hpp"
#include "refine_alignment.hpp"
#include "global.hpp"

using namespace bloomrepeats;

BOOST_AUTO_TEST_CASE (similar_aligner_n) {
    Strings seqs((3));
    seqs[0] = "ATTT";
    seqs[1] = "ANTT";
    seqs[2] = "ATTT";
    SimilarAligner::similar_aligner(seqs);
    BOOST_CHECK(seqs[0] == "ATTT");
    BOOST_CHECK(seqs[1] == "ANTT");
    BOOST_CHECK(seqs[2] == "ATTT");
}

BOOST_AUTO_TEST_CASE (similar_aligner_gap) {
    Strings seqs((3));
    seqs[0] = "ATGC";
    seqs[1] = "AGC";
    seqs[2] = "ATGC";
    SimilarAligner::similar_aligner(seqs, 1, 2, 5);
    BOOST_CHECK(seqs[0] == "ATGC");
    BOOST_CHECK(seqs[1] == "A-GC");
    BOOST_CHECK(seqs[2] == "ATGC");
}

BOOST_AUTO_TEST_CASE (similar_aligner_long_gap) {
    Strings seqs((3));
    seqs[0] = "ACCAGCTTTCGACCGCGGTGGCGATCGCGATATTAG";
    seqs[1] = "ACCAGCTGGTGGCGATCGCGATATTAG";
    seqs[2] = "ACCAGCTTTCGACCGCGGTGGCGATCGCGATATTAG";
    SimilarAligner::similar_aligner(seqs, 1, 2, 5);
    BOOST_CHECK(seqs[0] == "ACCAGCTTTCGACCGCGGTGGCGATCGCGATATTAG");
    BOOST_CHECK(seqs[1] == "ACCAGCT---------GGTGGCGATCGCGATATTAG");
    BOOST_CHECK(seqs[2] == "ACCAGCTTTCGACCGCGGTGGCGATCGCGATATTAG");
}

BOOST_AUTO_TEST_CASE (similar_aligner_empty) {
    Strings seqs((3));
    seqs[0] = "ATG";
    seqs[1] = "AG";
    SimilarAligner::similar_aligner(seqs, 1, 2, 5);
    BOOST_CHECK(seqs[0].size() == seqs[1].size());
    BOOST_CHECK(seqs[0].size() == seqs[2].size());
}

BOOST_AUTO_TEST_CASE (similar_aligner_end) {
    Strings seqs((2));
    seqs[0] = "ATG";
    seqs[1] = "AG";
    SimilarAligner::similar_aligner(seqs, 1, 2, 5);
    BOOST_CHECK(seqs[0] == "ATG");
    BOOST_CHECK(seqs[1] == "A-G");
}

BOOST_AUTO_TEST_CASE (similar_aligner_gap_repeat) {
    Strings seqs((2));
    seqs[0] = "CGAAT";
    seqs[1] = "CAAAT";
    SimilarAligner::similar_aligner(seqs, 1, 2, 5);
    BOOST_CHECK(seqs[0] == "CGAAT");
    BOOST_CHECK(seqs[1] == "CAAAT");
}

BOOST_AUTO_TEST_CASE (similar_aligner_refinement) {
    Strings seqs((2));
    seqs[0] = "CGGGCCGG";
    seqs[1] = "CGGCGCGG";
    SimilarAligner::similar_aligner(seqs, 1, 2, 5);
    refine_alignment(seqs);
    BOOST_CHECK(seqs[0] == "CGGGC-CGG");
    BOOST_CHECK(seqs[1] == "CGG-CGCGG");
}

BOOST_AUTO_TEST_CASE (similar_aligner_refinement_2) {
    Strings seqs((4));
    seqs[0] = "CGGGCCGG";
    seqs[1] = "CGGGCCGG";
    seqs[2] = "CGGCGCGG";
    seqs[3] = "CGGCGCGG";
    SimilarAligner::similar_aligner(seqs, 1, 2, 5);
    refine_alignment(seqs);
    BOOST_CHECK(seqs[0] == "CGGGC-CGG");
    BOOST_CHECK(seqs[1] == "CGGGC-CGG");
    BOOST_CHECK(seqs[2] == "CGG-CGCGG");
    BOOST_CHECK(seqs[3] == "CGG-CGCGG");
}

BOOST_AUTO_TEST_CASE (refinement_3) {
    Strings seqs((4));
    seqs[0] = "CCGG";
    seqs[1] = "CG-G";
    seqs[2] = "CG-G";
    seqs[3] = "CG-G";
    refine_alignment(seqs);
    BOOST_CHECK(seqs[0] == "CCGG");
    BOOST_CHECK(seqs[1] == "C-GG");
    BOOST_CHECK(seqs[2] == "C-GG");
    BOOST_CHECK(seqs[3] == "C-GG");
}

