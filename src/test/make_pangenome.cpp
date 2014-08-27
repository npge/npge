/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "MakePangenome.hpp"
#include "hit.hpp"
#include "CheckNoOverlaps.hpp"
#include "cast.hpp"

BOOST_AUTO_TEST_CASE (make_pangenome_TACG) {
    using namespace npge;
    std::string TACG;
    while (TACG.length() < 518) {
        TACG += "TACG";
    }
    SequencePtr s1 = boost::make_shared<InMemorySequence>(TACG);
    MakePangenome mp;
    mp.set_options("--out-stats=:null");
    mp.set_options("--min-fragment=100");
    mp.block_set()->add_sequence(s1);
    mp.run();
    BOOST_FOREACH (Block* block, *mp.block_set()) {
        BOOST_CHECK(!has_self_overlaps(block));
    }
    CheckNoOverlaps cno;
    cno.set_block_set(mp.block_set());
    cno.run();
    // check blocks found, not only unique part
    BOOST_CHECK(mp.block_set()->size() >= 2);
}

