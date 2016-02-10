/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Block.hpp"
#include "hit.hpp"

BOOST_AUTO_TEST_CASE (hit_main) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block block;
    block.insert(new Fragment(s1, 0, 5, 1));
    block.insert(new Fragment(s1, 5, 10, 1));
    BOOST_CHECK(has_self_overlaps(&block));
    BOOST_CHECK(block.alignment_length() == 6);
    fix_self_overlaps(&block);
    BOOST_CHECK(block.alignment_length() == 5);
}

BOOST_AUTO_TEST_CASE (hit_self_overlaps_rotation) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<CompactSequence>(std::string(200, 'A'));
    SequencePtr s2 = boost::make_shared<CompactSequence>(std::string(200, 'A'));
    Block block;
    block.insert(new Fragment(s1, 1, 10, 1));
    block.insert(new Fragment(s2, 5, 15, 1));
    block.insert(new Fragment(s1, 6, 15, 1));
    block.insert(new Fragment(s2, 10, 20, 1));
    BOOST_CHECK(has_self_overlaps(&block));
}

BOOST_AUTO_TEST_CASE (hit_main2) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block block;
    block.insert(new Fragment(s1, 0, 5, 1));
    block.insert(new Fragment(s1, 0, 5, 1));
    BOOST_CHECK(has_self_overlaps(&block));
    fix_self_overlaps(&block);
    BOOST_CHECK(block.empty());
}

BOOST_AUTO_TEST_CASE (hit_main3) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Block block;
    block.insert(new Fragment(s1, 0, 5, 1));
    block.insert(new Fragment(s1, 0, 5, -1));
    BOOST_CHECK(has_self_overlaps(&block));
    fix_self_overlaps(&block);
    BOOST_CHECK(block.size() == 2);
    BOOST_CHECK(block.alignment_length() == 3); // [0, 1, 2], [5, 4, 3]
}

