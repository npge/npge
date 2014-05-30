/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "SimilarAligner.hpp"
#include "refine_alignment.hpp"
#include "global.hpp"

using namespace npge;

BOOST_AUTO_TEST_CASE (similar_aligner_n) {
    Strings seqs((3));
    seqs[0] = "ATTT";
    seqs[1] = "ANTT";
    seqs[2] = "ATTT";
    SimilarAligner().similar_aligner(seqs);
    BOOST_CHECK(seqs[0] == "ATTT");
    BOOST_CHECK(seqs[1] == "ANTT");
    BOOST_CHECK(seqs[2] == "ATTT");
}

BOOST_AUTO_TEST_CASE (similar_aligner_gap) {
    Strings seqs((3));
    seqs[0] = "ATGC";
    seqs[1] = "AGC";
    seqs[2] = "ATGC";
    SimilarAligner().similar_aligner(seqs);
    BOOST_CHECK(seqs[0] == "ATGC");
    BOOST_CHECK(seqs[1] == "A-GC");
    BOOST_CHECK(seqs[2] == "ATGC");
}

BOOST_AUTO_TEST_CASE (similar_aligner_gap_2) {
    Strings seqs((2));
    seqs[0] = "CCCATATGG";
    seqs[1] = "CCATATCG";
    SimilarAligner().similar_aligner(seqs);
    refine_alignment(seqs);
    BOOST_CHECK(seqs[0] == "CCCATATGG");
    BOOST_CHECK(seqs[1] == "CC-ATATCG" || seqs[1] == "-CCATATCG" ||
                seqs[1] == "C-CATATCG");
}

BOOST_AUTO_TEST_CASE (similar_aligner_long_gap) {
    Strings seqs((3));
    seqs[0] = "ACCAGCTTTCGACCGCGGTGGCGATCGCGATATTAG";
    seqs[1] = "ACCAGCTGGTGGCGATCGCGATATTAG";
    seqs[2] = "ACCAGCTTTCGACCGCGGTGGCGATCGCGATATTAG";
    SimilarAligner().similar_aligner(seqs);
    BOOST_CHECK(seqs[0] == "ACCAGCTTTCGACCGCGGTGGCGATCGCGATATTAG");
    BOOST_CHECK(seqs[1] == "ACCAGCT---------GGTGGCGATCGCGATATTAG");
    BOOST_CHECK(seqs[2] == "ACCAGCTTTCGACCGCGGTGGCGATCGCGATATTAG");
}

BOOST_AUTO_TEST_CASE (similar_aligner_empty) {
    Strings seqs((3));
    seqs[0] = "ATG";
    seqs[1] = "AG";
    SimilarAligner().similar_aligner(seqs);
    BOOST_CHECK(seqs[0].size() == seqs[1].size());
    BOOST_CHECK(seqs[0].size() == seqs[2].size());
}

BOOST_AUTO_TEST_CASE (similar_aligner_end) {
    Strings seqs((2));
    seqs[0] = "ATG";
    seqs[1] = "AG";
    SimilarAligner().similar_aligner(seqs);
    BOOST_CHECK(seqs[0] == "ATG");
    BOOST_CHECK(seqs[1] == "A-G");
}

BOOST_AUTO_TEST_CASE (similar_aligner_gap_repeat) {
    Strings seqs((2));
    seqs[0] = "CGAAT";
    seqs[1] = "CAAAT";
    SimilarAligner().similar_aligner(seqs);
    BOOST_CHECK(seqs[0] == "CGAAT");
    BOOST_CHECK(seqs[1] == "CAAAT");
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

BOOST_AUTO_TEST_CASE (refinement_4) {
    Strings seqs((2));
    seqs[0] = "CCGGCC";
    seqs[1] = "CGG--C";
    refine_alignment(seqs);
    BOOST_CHECK(seqs[0] == "CCGGCC");
    BOOST_CHECK(seqs[1] == "C-GG-C");
}

BOOST_AUTO_TEST_CASE (refinement_5) {
    Strings seqs((3));
    seqs[0] = "-CCCCCC";
    seqs[1] = "CCCACCC";
    seqs[2] = "CCCTCCC";
    refine_alignment(seqs);
    BOOST_CHECK(seqs[0] == "CCC-CCC");
    BOOST_CHECK(seqs[1] == "CCCACCC");
    BOOST_CHECK(seqs[2] == "CCCTCCC");
}

BOOST_AUTO_TEST_CASE (similar_aligner_end_gap) {
    Strings seqs((2));
    seqs[0] = "GTTT";
    seqs[1] = "GTTTT";
    SimilarAligner().similar_aligner(seqs);
    BOOST_CHECK(seqs[0].length() == 5);
    BOOST_CHECK(seqs[0][4] != '-');
}

BOOST_AUTO_TEST_CASE (similar_aligner_bad) {
    Strings seqs((4));
    seqs[0] = "GCTATAAAGCAGCCTTCTTAGCTCACC";
    seqs[1] = "ACTTGATGTGCGGCTCGGGATATTTCA";
    seqs[2] = "CCCTCTCTGGGCAGGGCGAACATTAAA";
    seqs[3] = "TTGTAATGCTATTCCATAGTGAGATGA";
    SimilarAligner().similar_aligner(seqs);
}

BOOST_AUTO_TEST_CASE (similar_aligner_repeat_with_mismatch) {
    Strings seqs((2));
    seqs[0] = "AGAGCGGTTCCGGCGATTCCGTT";
    seqs[1] = "AGAGCGATTCCGTT";
    SimilarAligner().similar_aligner(seqs);
    // AGAGCGGTTCCGGCGATTCCGTT
    // AGAGCG******---ATTCCGTT
    // AGAGC-******--GATTCCGTT
    // AGAG--******-CGATTCCGTT
    // AGA---******GCGATTCCGTT
    // 01234567890123456
    BOOST_CHECK(seqs[0] == "AGAGCGGTTCCGGCGATTCCGTT");
    BOOST_CHECK(seqs[1].substr(6, 6) == "------");
}

