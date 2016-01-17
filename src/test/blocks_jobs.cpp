/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include <luabind/luabind.hpp>

#include "BlocksJobs.hpp"
#include "Meta.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

namespace npge {

class TestBlocksJobs : public BlocksJobs {
public:
    void process_block_impl(Block* b, ThreadData*) const {
        lua_State* L = meta()->L();
        luabind::object block(L, b);
        block["clear"](block);
    }
};

}

BOOST_AUTO_TEST_CASE (BlocksJobs_L) {
    using namespace npge;
    TestBlocksJobs tbj;
    SequencePtr seq(new InMemorySequence("TGAGATGCGGGCC"));
    tbj.block_set()->add_sequence(seq);
    for (int i = 0; i < 100; i++) {
        Block* b = new Block;
        b->insert(new Fragment(seq, 1, 2));
        tbj.block_set()->insert(b);
    }
    tbj.set_workers(10);
    tbj.run();
    BOOST_FOREACH (Block* b, *tbj.block_set()) {
        BOOST_CHECK(b->empty());
    }
}

