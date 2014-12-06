/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <boost/test/unit_test.hpp>

#include <lua.hpp>

#include "global.hpp"
#include "Meta.hpp"
#include "algo_lua.hpp"
#include "Decimal.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "FragmentCollection.hpp"
#include "AlignmentRow.hpp"
#include "Union.hpp"

BOOST_AUTO_TEST_CASE (MapAny_main) {
    using namespace npge;
    MapAny m1;
    Meta* meta = Meta::instance();
    lua_State* L = meta->L();
    //
    m1["bool"] = true;
    m1["int"] = 42;
    m1["Decimal"] = D(123.5);
    m1["string"] = std::string("string");
    Strings strings;
    strings.push_back("strings");
    // Strings must not be empty to be distinguishable
    // from other tables on Lua (typeless) side
    m1["Strings"] = strings;
    SequencePtr s(new InMemorySequence("ATGC"));
    Fragment* f = new Fragment(s, 1, 2);
    AlignmentRow* row = new CompactAlignmentRow;
    f->set_row(row);
    row->grow("T-G");
    Block* block = new Block;
    block->insert(f);
    BlockSetPtr bs = new_bs();
    bs->insert(block);
    m1["Sequence"] = s;
    m1["Fragment"] = f;
    m1["AlignmentRow"] = row;
    m1["Block"] = block;
    m1["BlockSet"] = bs;
    Blocks blocks;
    blocks.push_back(block);
    Fragments fragments;
    fragments.push_back(f);
    m1["Blocks"] = blocks;
    m1["Fragments"] = fragments;
    typedef boost::shared_ptr<SetFc> SetFcPtr;
    SetFcPtr set_fc(new SetFc);
    set_fc->add_bs(*bs);
    set_fc->prepare();
    typedef boost::shared_ptr<VectorFc> VectorFcPtr;
    VectorFcPtr vector_fc(new VectorFc);
    vector_fc->add_bs(*bs);
    vector_fc->prepare();
    m1["SetFc"] = set_fc;
    m1["VectorFc"] = vector_fc;
    SharedProcessor u = meta->get("Union");
    m1["Processor"] = u.get();
    //
    typedef luabind::default_converter<npge::MapAny> dcM;
    dcM dcm;
    dcm.to(L, m1);
    BOOST_CHECK(dcm.compute_score(L, -1) >= 0);
    MapAny m2 = dcm.from(L, -1);
    lua_pop(L, 1);
    //
    BOOST_CHECK(m2.size() == m1.size());
    BOOST_CHECK(m2["bool"].as<bool>() == true);
    BOOST_CHECK(m2["int"].as<int>() == 42);
    BOOST_CHECK(m2["Decimal"].as<Decimal>() == D(123.5));
    BOOST_CHECK(m2["string"].as<std::string>() == "string");
    BOOST_CHECK(m2["Strings"].as<Strings>()[0] == "strings");
    BOOST_CHECK(m2["Sequence"].as<SequencePtr>() == s);
    BOOST_CHECK(m2["Fragment"].as<Fragment*>() == f);
    BOOST_CHECK(m2["AlignmentRow"].as<AlignmentRow*>() == row);
    BOOST_CHECK(m2["Block"].as<Block*>() == block);
    BOOST_CHECK(m2["BlockSet"].as<BlockSetPtr>() == bs);
    BOOST_CHECK(m2["Processor"].as<Processor*>() == u.get());
    BOOST_CHECK(m2["Blocks"].as<Blocks>()[0] == block);
    BOOST_CHECK(m2["Fragments"].as<Fragments>()[0] == f);
    BOOST_CHECK(m2["SetFc"].as<SetFcPtr>() == set_fc);
    BOOST_CHECK(m2["VectorFc"].as<VectorFcPtr>() == vector_fc);
}

