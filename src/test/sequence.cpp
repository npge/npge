/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

BOOST_AUTO_TEST_CASE (Sequence_main) {
    using namespace npge;
    InMemorySequence seq("TGGTCCGAGATGCGGGCCCGTAAGCTTACATACAGG");
    BOOST_REQUIRE(seq.size() == 36);
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGG");
    BOOST_REQUIRE(s1->size() == 3);
    size_t fragments_number = 0;
    Fragment f(s1);
    s1->make_first_fragment(f, 2);
    while (s1->next_fragment(f)) {
        fragments_number += 1;
    }
    BOOST_CHECK(fragments_number == 4);
}

BOOST_AUTO_TEST_CASE (Sequence_n) {
    using namespace npge;
    std::string seq_str = "TGGTCNGAGATGCGG";
    InMemorySequence in_mem_seq(seq_str);
    BOOST_CHECK(in_mem_seq.size() == seq_str.size());
    BOOST_CHECK(in_mem_seq.contents() == seq_str);
    CompactSequence compact_seq(seq_str);
    BOOST_CHECK(compact_seq.size() == seq_str.size());
    BOOST_CHECK(compact_seq.contents() == seq_str);
}

BOOST_AUTO_TEST_CASE (Sequence_first_ori) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGG");
    Fragment f(s1);
    s1->make_first_fragment(f, 2);
    while (s1->next_fragment(f)) {
        BOOST_CHECK(f.ori() == Sequence::FIRST_ORI);
        break;
    }
}

BOOST_AUTO_TEST_CASE (Sequence_filtering) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>(" ---ATG--caggacg..");
    BOOST_REQUIRE(s1->size() == 10);
    Fragment f(s1, 0, 9);
    BOOST_CHECK(f.str() == "ATGCAGGACG");
    BOOST_CHECK(f.at(-1) == 'G');
    BOOST_CHECK(s1->contents() == "ATGCAGGACG");
}

BOOST_AUTO_TEST_CASE (Sequence_compact_main) {
    using namespace npge;
    CompactSequence seq("TGGTCCGAGATGCGGGCCCGTAAGCTTACATACAGG");
    BOOST_REQUIRE(seq.size() == 36);
    SequencePtr s1 = boost::make_shared<CompactSequence>("TGG");
    BOOST_REQUIRE(s1->size() == 3);
    size_t fragments_number = 0;
    Fragment f(s1);
    s1->make_first_fragment(f, 2);
    while (s1->next_fragment(f)) {
        fragments_number += 1;
    }
    BOOST_CHECK(fragments_number == 4);
}

BOOST_AUTO_TEST_CASE (Sequence_compact_filtering) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<CompactSequence>(" ---ATG--caggacg..");
    BOOST_REQUIRE(s1->size() == 10);
    Fragment f(s1, 0, 9);
    BOOST_CHECK(f.str() == "ATGCAGGACG");
    BOOST_CHECK(f.at(-1) == 'G');
}

BOOST_AUTO_TEST_CASE (Sequence_genome_chromosome) {
    using namespace npge;
    CompactSequence s("ATG");
    s.set_name("123");
    BOOST_REQUIRE(s.name() == "123");
    BOOST_CHECK(s.genome() == "");
    BOOST_CHECK(s.chromosome() == "");
    s.set_name("abc&chr1");
    BOOST_CHECK(s.genome() == "");
    BOOST_CHECK(s.chromosome() == "");
    s.set_name("abc&chr1&c");
    BOOST_CHECK(s.genome() == "abc");
    BOOST_CHECK(s.chromosome() == "chr1");
    s.set_name("abc&chr1&zzz");
    BOOST_CHECK(s.genome() == "");
    BOOST_CHECK(s.chromosome() == "");
    s.set_name("abc&chr1&circular");
    BOOST_CHECK(s.genome() == "");
    BOOST_CHECK(s.chromosome() == "");
    s.set_name("abc&chr1&linear");
    BOOST_CHECK(s.genome() == "");
    BOOST_CHECK(s.chromosome() == "");
    bool thrown = false;
    try {
        s.set_name("abc&chr1&linear");
        bool circular = s.circular();
    } catch (...) {
        thrown = true;
    }
    BOOST_CHECK(thrown);
    thrown = false;
    try {
        s.set_name("abc&chr1&c");
        bool circular = s.circular();
        BOOST_CHECK(circular);
    } catch (...) {
        thrown = true;
    }
    BOOST_CHECK(!thrown);
}

BOOST_AUTO_TEST_CASE (Sequence_consensus_of_block) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<CompactSequence>("CAGGACGG");
    SequencePtr s2 = boost::make_shared<CompactSequence>("CAGGAAG-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("CTGGACG-");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    Block block;
    block.insert(f1);
    block.insert(f2);
    block.insert(f3);
    block.set_name("myblock");
    BOOST_WARN(block.consensus_string() == "CAGGACGG");
    InMemorySequence consensus("");
    consensus.set_name("myname");
    consensus.set_description("mydescr");
    consensus.set_block(&block);
    BOOST_CHECK(consensus.block() == &block);
    BOOST_CHECK(consensus.contents() == block.consensus_string());
    BOOST_CHECK(consensus.name() == "myname");
    BOOST_CHECK(consensus.description() == "mydescr");
}

BOOST_AUTO_TEST_CASE (Sequence_consensus_of_block_empty_name) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<CompactSequence>("CAGGACGG");
    SequencePtr s2 = boost::make_shared<CompactSequence>("CAGGAAG-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("CTGGACG-");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    Block block;
    block.insert(f1);
    block.insert(f2);
    block.insert(f3);
    block.set_name("myblock");
    BOOST_WARN(block.consensus_string() == "CAGGACGG");
    InMemorySequence consensus("");
    consensus.set_name("");
    consensus.set_description("");
    consensus.set_block(&block);
    BOOST_CHECK(consensus.block() == &block);
    BOOST_CHECK(consensus.contents() == block.consensus_string());
    BOOST_CHECK(consensus.name() == "myblock");
}

