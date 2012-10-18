/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>
#include <boost/foreach.hpp>

#include "AddBlocks.hpp"
#include "BlockSet.hpp"
#include "Alignment.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

AddBlocks::AddBlocks(bool keep_alignment):
    keep_alignment_(keep_alignment)
{ }

void AddBlocks::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("in-blocks", po::value<Files>()->required(),
     "input fasta file(s) with blocks")
    ("keep-alignment", po::bool_switch(),
     "import alignment (not only start and stop positions)")
    ("no-keep-alignment", po::bool_switch(), "do not import alignment")
   ;
}

void AddBlocks::apply_options_impl(const po::variables_map& vm) {
    set_files(vm["in-blocks"].as<Files>());
    if (vm["keep-alignment"].as<bool>() && vm["no-keep-alignment"].as<bool>()) {
        throw Exception("both 'keep-alignment' and "
                        "'no-keep-alignment' specified");
    }
    if (vm["keep-alignment"].as<bool>()) {
        set_keep_alignment(true);
    }
    if (vm["no-keep-alignment"].as<bool>()) {
        set_keep_alignment(false);
    }
}

bool AddBlocks::run_impl() const {
    int size_before = block_set()->size();
    BOOST_FOREACH (std::string file_name, files()) {
        std::ifstream input_file(file_name.c_str());
        BlockSetFastaReader reader(*block_set(), input_file, keep_alignment());
        reader.read_all_sequences();
    }
    return block_set()->size() > size_before;
}

const char* AddBlocks::name_impl() const {
    return "Input block set";
}

}

