/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/assert.hpp>

#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "AddSequences.hpp"
#include "AddBlocks.hpp"
#include "CleanUp.hpp"
#include "UniqueNames.hpp"
#include "Output.hpp"
#ifndef NDEBUG
#include "Connector.hpp"
#include "OverlapsResolver.hpp"
#endif
#include "po.hpp"
#include "Exception.hpp"

using namespace bloomrepeats;

int main(int argc, char** argv) {
    po::options_description desc("Options");
    add_general_options(desc);
    AddSequences adder;
    adder.add_options(desc);
    AddBlocks blocks_adder;
    blocks_adder.add_options(desc);
    po::positional_options_description pod;
    Output output;
    output.add_options(desc);
    CleanUp cleanup;
    cleanup.add_options(desc);
    po::variables_map vm;
    int error = read_options(argc, argv, vm, desc, pod);
    if (error) {
        return error;
    }
    try {
        adder.apply_options(vm);
        blocks_adder.apply_options(vm);
        cleanup.apply_options(vm);
        output.apply_options(vm);
    } catch (Exception& e) {
        std::cerr << argv[0] << ": " << e.what() << std::endl;
        return 255;
    }
    BlockSetPtr pangenome = boost::make_shared<BlockSet>();
    adder.apply(pangenome);
    blocks_adder.apply(pangenome);
    cleanup.apply(pangenome);
#ifndef NDEBUG
    Connector connector;
    connector.apply(pangenome);
    OverlapsResolver resolver;
    resolver.set_block_set(pangenome);
    BOOST_ASSERT(!resolver.overlaps());
#endif
    UniqueNames names;
    names.apply(pangenome);
    output.apply(pangenome);
}

