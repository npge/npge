/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>

#include "Pipe.hpp"
#include "Connector.hpp"
#include "FragmentsExpander.hpp"
#include "BlocksExpander.hpp"
#include "Joiner.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

BOOST_AUTO_TEST_CASE (Pipe_main) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggCC");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("tGGtccgagcggacggcc");
    Block* b1 = new Block();
    b1->insert(new Fragment(s1, 1, 2));
    b1->insert(new Fragment(s2, 1, 2));
    Block* b2 = new Block();
    b2->insert(new Fragment(s1, 16, 17));
    BlockSetPtr block_set = new_bs();
    block_set->insert(b1);
    block_set->insert(b2);
    Pipe expander;
    expander.add(new Connector);
    expander.add(new FragmentsExpander).add(new BlocksExpander).add(new Joiner);
    expander.apply(block_set);
    BOOST_REQUIRE(block_set->size() == 1);
    BOOST_CHECK(block_set->front()->size() == 2);
    BOOST_CHECK(block_set->front()->front()->length() == s1->size());
}

