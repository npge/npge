/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "OriByMajority.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

BOOST_AUTO_TEST_CASE (OriByMajority_main) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGGTCCGAGCGGACGGCC");
    Fragment* f1 = new Fragment(s1, 2, 6, 1);
    Fragment* f2 = new Fragment(s1, 9, 13, -1);
    Block block;
    block.insert(f1);
    block.insert(f2);
    OriByMajority obm;
    obm.apply_to_block(&block);
    BOOST_CHECK(f1->ori() == 1);
    block.insert(new Fragment(s1, 9, 13, -1));
    obm.apply_to_block(&block);
    BOOST_CHECK(f1->ori() == -1);
}

