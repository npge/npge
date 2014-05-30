/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>

#include "Sequence.hpp"
#include "AlignmentRow.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "Filter.hpp"
#include "SizeLimits.hpp"

BOOST_AUTO_TEST_CASE (Filter_good_block) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("AAA");
    Filter filter;
    allow_everything(&filter);
    filter.set_options("--min-fragment=100");
    filter.set_options("--min-block=2");
    boost::scoped_ptr<Block> block(new Block);
    block->insert(new Fragment(s1, 0, 0));
    block->insert(new Fragment(s1, 1, 1));
    BOOST_CHECK(!filter.is_good_block(block.get()));
    filter.set_options("--min-fragment=1");
    BOOST_CHECK(filter.is_good_block(block.get()));
}

BOOST_AUTO_TEST_CASE (Filter_good_blocks) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TAGTCCG-");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGTT-CGT");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TG---CG-");
    Block b;
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    new MapAlignmentRow("TAGTCCG-", f1);
    b.insert(f1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TGTT-CGT", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    new MapAlignmentRow("TG---CG-", f3);
    b.insert(f3);
    Filter filter;
    allow_everything(&filter);
    filter.set_opt_value("min-fragment", 2);
    filter.set_opt_value("min-block", 2);
    filter.set_opt_value("min-identity", 0.99);
    filter.set_opt_value("max-gaps", 0.01);
    std::vector<Block*> good_blocks;
    filter.find_good_subblocks(&b, good_blocks);
    BOOST_CHECK(good_blocks.size() == 1);
    if (good_blocks.size() == 1) {
        Block* gb = good_blocks[0];
        BOOST_CHECK(gb->size() == 3);
        BOOST_CHECK(gb->alignment_length() == 2);
        BOOST_CHECK(gb->front()->str() == "CG");
    }
    BOOST_FOREACH (Block* block, good_blocks) {
        delete block;
    }
}

BOOST_AUTO_TEST_CASE (Filter_good_blocks2) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TATTCCG-");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGTT-CGT");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TGT--CG-");
    Block b;
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    new MapAlignmentRow("TATTCCG-", f1);
    b.insert(f1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TGTT-CGT", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    new MapAlignmentRow("TGT--CG-", f3);
    b.insert(f3);
    Filter filter;
    allow_everything(&filter);
    filter.set_opt_value("min-fragment", 3);
    filter.set_opt_value("min-block", 2);
    filter.set_opt_value("min-identity", 0.6);
    filter.set_opt_value("max-gaps", 0.01);
    std::vector<Block*> good_blocks;
    filter.find_good_subblocks(&b, good_blocks);
    BOOST_CHECK(good_blocks.size() == 1);
    BOOST_CHECK(good_blocks.size() == 1);
    if (good_blocks.size() == 1) {
        Block* gb = good_blocks[0];
        BOOST_CHECK(gb->size() == 3);
        BOOST_CHECK(gb->alignment_length() == 3);
        BOOST_CHECK(gb->front()->str() == "TAT" ||
                    gb->front()->str() == "TGT");
    }
    BOOST_FOREACH (Block* block, good_blocks) {
        delete block;
    }
}

BOOST_AUTO_TEST_CASE (Filter_good_blocks3) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TATTCCG-");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TGTTACGT");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("TGT--CG-");
    Block b;
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    new MapAlignmentRow("TATTCCG-", f1);
    b.insert(f1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TGTTACGT", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    new MapAlignmentRow("TGT--CG-", f3);
    b.insert(f3);
    Filter filter;
    allow_everything(&filter);
    filter.set_opt_value("min-fragment", 1);
    filter.set_opt_value("min-block", 2);
    filter.set_opt_value("min-identity", 0.99);
    filter.set_opt_value("max-gaps", 0.51);
    std::vector<Block*> good_blocks;
    filter.find_good_subblocks(&b, good_blocks);
    BOOST_CHECK(good_blocks.size() == 3);
    BOOST_FOREACH (Block* block, good_blocks) {
        delete block;
    }
}

BOOST_AUTO_TEST_CASE (Filter_good_blocks_expand) {
    using namespace npge;
    SequencePtr s1 = boost::make_shared<InMemorySequence>("TGTTCCG");
    SequencePtr s2 = boost::make_shared<InMemorySequence>("TATTCC-");
    SequencePtr s3 = boost::make_shared<InMemorySequence>("-ATTCCG");
    Block b;
    Fragment* f1 = new Fragment(s1, 0, s1->size() - 1);
    new MapAlignmentRow("TGTTCCG", f1);
    b.insert(f1);
    Fragment* f2 = new Fragment(s2, 0, s2->size() - 1);
    new MapAlignmentRow("TATTCC-", f2);
    b.insert(f2);
    Fragment* f3 = new Fragment(s3, 0, s3->size() - 1);
    new MapAlignmentRow("-ATTCCG", f3);
    b.insert(f3);
    Filter filter;
    allow_everything(&filter);
    filter.set_opt_value("min-fragment", 1);
    filter.set_opt_value("min-block", 2);
    filter.set_opt_value("min-identity", 0.99);
    filter.set_opt_value("max-gaps", 0.01);
    std::vector<Block*> good_blocks;
    filter.find_good_subblocks(&b, good_blocks);
    BOOST_CHECK(good_blocks.size() == 1);
    if (good_blocks.size() == 1) {
        Block* gb = good_blocks[0];
        BOOST_CHECK(gb->size() == 3);
        BOOST_CHECK(gb->alignment_length() == 4);
        BOOST_CHECK(gb->front()->str() == "TTCC");
    }
    BOOST_FOREACH (Block* block, good_blocks) {
        delete block;
    }
}

