/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "CleanUp.hpp"
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
    Sequence::add_input_options(desc);
    po::positional_options_description pod;
    BlockSet::add_output_options(desc);
    desc.add_options()
    ("pangenome", po::value<std::string>()->required(),
     "input file with existing pangenome");
    CleanUp cleanup;
    cleanup.add_options(desc);
    po::variables_map vm;
    int error = read_options(argc, argv, vm, desc, pod);
    if (error) {
        return error;
    }
    try {
        cleanup.apply_options(vm);
    } catch (Exception& e) {
        std::cerr << argv[0] << ": " << e.what() << std::endl;
        return 255;
    }
    std::vector<SequencePtr> seqs;
    Sequence::read_all_files(vm, seqs);
    BlockSetPtr pangenome = boost::make_shared<BlockSet>();
    pangenome->add_sequences(seqs);
    std::ifstream pangenome_file(vm["pangenome"].as<std::string>().c_str());
    pangenome_file >> *pangenome;
    cleanup.apply(pangenome);
#ifndef NDEBUG
    Connector connector;
    connector.apply(pangenome);
    OverlapsResolver resolver;
    resolver.set_block_set(pangenome);
    BOOST_ASSERT(!resolver.overlaps());
#endif
    pangenome->set_unique_block_names();
    pangenome->make_output(vm);
}

