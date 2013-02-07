/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Block.hpp"
#include "Fragment.hpp"

BOOST_AUTO_TEST_CASE (Sequence_main) {
    using namespace bloomrepeats;
    InMemorySequence seq("tggtccgagatgcgggcccgtaagcttacatacagg");
    BOOST_REQUIRE(seq.size() == 36);
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgg");
    BOOST_REQUIRE(s1->size() == 3);
    size_t fragments_number = 0;
    Fragment f(s1);
    s1->make_first_fragment(f, 2);
    while (s1->next_fragment(f)) {
        fragments_number += 1;
    }
    BOOST_CHECK(fragments_number == 4);
}

BOOST_AUTO_TEST_CASE (Sequence_first_ori) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tgg");
    Fragment f(s1);
    s1->make_first_fragment(f, 2);
    while (s1->next_fragment(f)) {
        BOOST_CHECK(f.ori() == Sequence::FIRST_ORI);
        break;
    }
}

BOOST_AUTO_TEST_CASE (Sequence_filtering) {
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<InMemorySequence>(" ---ATG--caggacg..");
    BOOST_REQUIRE(s1->size() == 10);
    Fragment f(s1, 0, 9);
    BOOST_CHECK(f.str() == "atgcaggacg");
    BOOST_CHECK(f.at(-1) == 'g');
    BOOST_CHECK(s1->contents() == "atgcaggacg");
}

BOOST_AUTO_TEST_CASE (Sequence_compact_main) {
    using namespace bloomrepeats;
    CompactSequence seq("tggtccgagatgcgggcccgtaagcttacatacagg");
    BOOST_REQUIRE(seq.size() == 36);
    SequencePtr s1 = boost::make_shared<CompactSequence>("tgg");
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
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<CompactSequence>(" ---ATG--caggacg..");
    BOOST_REQUIRE(s1->size() == 10);
    Fragment f(s1, 0, 9);
    BOOST_CHECK(f.str() == "atgcaggacg");
    BOOST_CHECK(f.at(-1) == 'g');
}

BOOST_AUTO_TEST_CASE (Sequence_genome_chromosome) {
    using namespace bloomrepeats;
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
    using namespace bloomrepeats;
    SequencePtr s1 = boost::make_shared<CompactSequence>("caggacgg");
    SequencePtr s2 = boost::make_shared<CompactSequence>("caggaag-");
    SequencePtr s3 = boost::make_shared<CompactSequence>("ctggacg-");
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    Block block;
    block.insert(f1);
    block.insert(f2);
    block.insert(f3);
    BOOST_WARN(block.consensus_string() == "caggacgg");
    InMemorySequence consensus("");
    consensus.set_name("myname");
    consensus.set_description("mydescr");
    consensus.set_block(&block);
    BOOST_CHECK(consensus.block() == &block);
    BOOST_CHECK(consensus.contents() == block.consensus_string());
    BOOST_CHECK(consensus.name() == "myname");
    BOOST_CHECK(consensus.description() == "mydescr");
}

