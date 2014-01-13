/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <utility>
#include <boost/test/unit_test.hpp>

#include "PairAligner.hpp"
#include "Sequence.hpp"

BOOST_AUTO_TEST_CASE (PairAligner_main) {
    using namespace bloomrepeats;
    std::string s1("gaacaggcttgtTtatttttacgttccctctacgccgctccGaacgtgagactct");
    std::string s2("gaacaggcttgtAtatttttacgttccctctacgccgctccCaacgtgagactct");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(2);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s1.size() - 1);
    BOOST_CHECK(s2_last == s2.size() - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_n_positive) {
    using namespace bloomrepeats;
    std::string s1("gaacaggcttgtNtNttt");
    std::string s2("gaacaggcttgtAtNttt");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(2);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s1.size() - 1);
    BOOST_CHECK(s2_last == s2.size() - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_n_negative) {
    using namespace bloomrepeats;
    std::string s1("gaacaggcttgtNtNttt");
    std::string s2("gaacaggcttgtAtNttt");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last < s1.size() - 1);
    BOOST_CHECK(s2_last < s2.size() - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_many_substitutions) {
    using namespace bloomrepeats;
    std::string s1("gaacaggcttgtTtatttttacgAtccctctacgccgctccGa");
    std::string s2("gaacaggcttgtAtatttttacgTtccctctacgccgctccCa");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(2);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s1.size() - 2 - 1);
    BOOST_CHECK(s2_last == s2.size() - 2 - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_gaps) {
    using namespace bloomrepeats;
    std::string s1("gaacaggcttgt-tatgattacgatccctctacgccgctccGa");
    std::string s2("gaacaggcttgtatatgattacg-tccctctacgccgctccCa");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    BOOST_CHECK(s2 == "GAACAGGCTTGTATATGATTACGTCCCTCTACGCCGCTCCCA");
    PairAligner aligner(1, 5, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == 23 - 1 - 1);
    BOOST_CHECK(s2_last == 23 - 1);
    //
    aligner.set_gap_range(2);
    BOOST_REQUIRE(aligner.max_errors() == 2);
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s1.size() - 2 - 1);
    BOOST_CHECK(s2_last == s2.size() - 2 - 1);
    //
    aligner.set_gap_range(1);
    BOOST_REQUIRE(aligner.max_errors() == 2);
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s1.size() - 2 - 1); // gaps on different seqs
    BOOST_CHECK(s2_last == s2.size() - 2 - 1);
    //
    aligner.set_gap_range(0);
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == 12 - 1);
    BOOST_CHECK(s2_last == 12 - 1);
    //
    aligner.set_gap_range(5);
    aligner.set_gap_penalty(2);
    aligner.set_max_errors(2);
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == 23 - 1 - 1);
    BOOST_CHECK(s2_last == 23 - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_gaps_gaps) {
    using namespace bloomrepeats;
    std::string s1("GAACAGGCTTGT--GTTAT");
    std::string s2("GAACAGGCTTGTAAGTTAT");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(2, 5, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s1.size() - 1);
    BOOST_CHECK(s2_last == s2.size() - 1);
    //
    aligner.set_gap_range(1);
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == 12 - 1);
    BOOST_CHECK(s2_last == 12 - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_test_3) {
    using namespace bloomrepeats;
    std::string s1("GAACAG-CTTGT--GTTAT");
    std::string s2("GA-CAGGCT-GTAAGTT-T");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(6, 3, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s1.size() - 1);
    BOOST_CHECK(s2_last == s2.size() - 1);
    //
    aligner.set_max_errors(1);
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == 6 - 1);
    BOOST_CHECK(s2_last == 6 - 1 - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_alignment) {
    using namespace bloomrepeats;
    std::string s1("GAACAG-CTTGT--GTTAT");
    std::string s2("GA-CAGGCT-GTAAGTT-T");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(6, 3, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    std::string s1_str, s2_str;
    std::vector<std::pair<int, int> > alignment;
    aligner.align(s1_last, s2_last, &s1_str, &s2_str, &alignment);
    BOOST_CHECK(s1_str == "GAACAG-CTTGT--GTTAT");
    BOOST_CHECK(s2_str == "GA-CAGGCT-GTAAGTT-T");
    BOOST_CHECK(alignment[0] == std::make_pair(0, 0));
    BOOST_CHECK(alignment[2] == std::make_pair(2, -1));
    BOOST_CHECK(alignment[3] == std::make_pair(3, 2));
    BOOST_CHECK(alignment.size() == 19);
}

BOOST_AUTO_TEST_CASE (PairAligner_very_short_true) {
    using namespace bloomrepeats;
    std::string s1("G");
    std::string s2("G");
    PairAligner aligner(0, 5);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    std::string s1_str, s2_str;
    std::vector<std::pair<int, int> > alignment;
    aligner.align(s1_last, s2_last, &s1_str, &s2_str, &alignment);
    BOOST_CHECK(s1_last == 0);
    BOOST_CHECK(s2_last == 0);
    BOOST_CHECK(s1_str == "G");
    BOOST_CHECK(s2_str == "G");
    BOOST_CHECK(alignment.size() == 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_very_short_false) {
    using namespace bloomrepeats;
    std::string s1("G");
    std::string s2("A");
    PairAligner aligner(0, 5);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    std::string s1_str, s2_str;
    std::vector<std::pair<int, int> > alignment;
    aligner.align(s1_last, s2_last, &s1_str, &s2_str, &alignment);
    BOOST_CHECK(s1_last == -1);
    BOOST_CHECK(s2_last == -1);
    BOOST_CHECK(s1_str == "");
    BOOST_CHECK(s2_str == "");
    BOOST_CHECK(alignment.size() == 0);
}

BOOST_AUTO_TEST_CASE (PairAligner_empty) {
    using namespace bloomrepeats;
    std::string s1("");
    std::string s2("");
    PairAligner aligner(0, 5);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    std::string s1_str, s2_str;
    std::vector<std::pair<int, int> > alignment;
    aligner.align(s1_last, s2_last, &s1_str, &s2_str, &alignment);
    BOOST_CHECK(s1_last == -1);
    BOOST_CHECK(s2_last == -1);
    BOOST_CHECK(s1_str == "");
    BOOST_CHECK(s2_str == "");
    BOOST_CHECK(alignment.size() == 0);
}

BOOST_AUTO_TEST_CASE (PairAligner_bad_alignment) {
    using namespace bloomrepeats;
    std::string s1("GAACAG-CTTGT--GTTAT");
    std::string s2("GA-CAGGCT-GTAAGTT-T");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(1, 1, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    std::string s1_str, s2_str;
    std::vector<std::pair<int, int> > alignment;
    aligner.align(s1_last, s2_last, &s1_str, &s2_str, &alignment);
    BOOST_CHECK(s1_str == "GAACAG");
    BOOST_CHECK(s2_str == "GA-CAG");
    BOOST_REQUIRE(alignment.size() == 6);
    BOOST_CHECK(alignment[0] == std::make_pair(0, 0));
    BOOST_CHECK(alignment[2] == std::make_pair(2, -1));
    BOOST_CHECK(alignment[3] == std::make_pair(3, 2));
}

BOOST_AUTO_TEST_CASE (PairAligner_alignment_custom_gap) {
    using namespace bloomrepeats;
    std::string s1("GAACAG-CTTGT--GTTAT");
    std::string s2("GA-CAGGCT-GTAAGTT-T");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(6, 3, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    std::string s1_str, s2_str;
    std::vector<std::pair<int, int> > alignment;
    aligner.align(s1_last, s2_last, &s1_str, &s2_str, &alignment, '.');
    BOOST_CHECK(s1_str == "GAACAG.CTTGT..GTTAT");
    BOOST_CHECK(s2_str == "GA.CAGGCT.GTAAGTT.T");
}

BOOST_AUTO_TEST_CASE (PairAligner_tail) {
    using namespace bloomrepeats;
    std::string s1("TTCCGGTGCTGCGaggga");
    std::string s2("TTCCGGTGCTGCGcctct");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(5, 5, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    {
        std::string s1_str, s2_str;
        aligner.set_no_tail(false);
        aligner.align(s1_last, s2_last, &s1_str, &s2_str);
        BOOST_CHECK(s1_str.size() == 18); // gaps or mismatches
        BOOST_CHECK(s2_str.size() == 18); // gaps or mismatches
    }
    {
        std::string s1_str, s2_str;
        aligner.set_no_tail(true);
        aligner.align(s1_last, s2_last, &s1_str, &s2_str);
        BOOST_CHECK(s1_str == "TTCCGGTGCTGCG");
        BOOST_CHECK(s2_str == "TTCCGGTGCTGCG");
    }
}

BOOST_AUTO_TEST_CASE (PairAligner_cols_less_than_rows) {
    using namespace bloomrepeats;
    std::string s1("gaacaggcttgtTtatttttacgttccctctacgccgctccGaacgtgagactct");
    std::string s2("gaacaggcttgtAtatttttacg");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(2, 5, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s2.size() - 1);
    BOOST_CHECK(s2_last == s2.size() - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_rows_less_than_cols) {
    using namespace bloomrepeats;
    std::string s2("gaacaggcttgtTtatttttacgttccctctacgccgctccGaacgtgagactct");
    std::string s1("gaacaggcttgtAtatttttacg");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(2, 5, 1);
    aligner.set_first(s1.c_str(), s1.size());
    aligner.set_second(s2.c_str(), s2.size());
    int s1_last, s2_last;
    aligner.align(s1_last, s2_last);
    BOOST_CHECK(s1_last == s1.size() - 1);
    BOOST_CHECK(s2_last == s1.size() - 1);
}

BOOST_AUTO_TEST_CASE (PairAligner_aligned) {
    using namespace bloomrepeats;
    PairAligner a(2, 2, 1);
    bool old_no_tail = a.no_tail();
    BOOST_CHECK(a.aligned("GAAC", "GAAC"));
    BOOST_CHECK(a.aligned("GAAC", "GA"));
    BOOST_CHECK(a.aligned("GAAC", "GTAC"));
    BOOST_CHECK(a.aligned("GAAC", "GTAT"));
    BOOST_CHECK(!a.aligned("GAAC", "GTTT"));
    BOOST_CHECK(!a.aligned("GAAC", "G"));
    BOOST_CHECK(a.aligned("GG", "GA"));
    BOOST_CHECK(a.aligned("GG", "AA"));
    a.set_max_errors(0);
    BOOST_CHECK(a.aligned("GA", "GA"));
    BOOST_CHECK(!a.aligned("GG", "GA"));
    BOOST_CHECK(a.aligned("G", "G"));
    BOOST_CHECK(a.aligned("A", "A"));
    BOOST_CHECK(a.aligned("T", "T"));
    BOOST_CHECK(a.aligned("C", "C"));
    BOOST_CHECK(!a.aligned("G", "A"));
    BOOST_CHECK(!a.aligned("T", "A"));
    BOOST_CHECK(a.aligned("", ""));
    BOOST_CHECK(a.no_tail() == old_no_tail);
}

BOOST_AUTO_TEST_CASE (PairAligner_aligned_last) {
    using namespace bloomrepeats;
    std::string s1("TTCCGGTGCTGCGaggga");
    std::string s2("TTCCGGTGCTGCGcctct");
    Sequence::to_atgcn(s1);
    Sequence::to_atgcn(s2);
    PairAligner aligner(5, 5, 1);
    int s1_last, s2_last;
    aligner.set_no_tail(false);
    BOOST_CHECK(aligner.aligned(s1, s2, &s1_last, &s2_last));
    BOOST_CHECK(s1_last == s1.size() - 1);
    BOOST_CHECK(s2.size() - 5 - 1 <= s2_last && s2_last <= s2.size() - 1);
    aligner.set_no_tail(true);
    BOOST_CHECK(aligner.aligned(s1, s2, &s1_last, &s2_last));
    BOOST_CHECK(s1_last == s1.size() - 5 - 1);
    BOOST_CHECK(s2_last == s2.size() - 5 - 1);
}

