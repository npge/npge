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
    keep_alignment_(keep_alignment),
    row_type_("compact")
{ }

void AddBlocks::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("in-blocks", po::value<Files>()->required(),
     "input fasta file(s) with blocks")
    ("keep-alignment", po::bool_switch(),
     "import alignment (not only start and stop positions)")
    ("no-keep-alignment", po::bool_switch(), "do not import alignment")
    ("row-type", po::value<std::string>()->default_value(row_type()),
     "way of storing alignments in memory ('map' or 'compact')");
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
    std::string row_type = vm["row-type"].as<std::string>();
    if (row_type != "map" && row_type != "compact") {
        throw Exception("'row-type' must be 'map' or 'compact'");
    }
    set_row_type(row_type);
}

bool AddBlocks::run_impl() const {
    int size_before = block_set()->size();
    BOOST_FOREACH (std::string file_name, files()) {
        std::ifstream input_file(file_name.c_str());
        RowType type = (row_type() == "map") ? MAP_ROW : COMPACT_ROW;
        BlockSetFastaReader reader(*block_set(), input_file,
                                   keep_alignment(), type);
        reader.read_all_sequences();
    }
    return block_set()->size() > size_before;
}

const char* AddBlocks::name_impl() const {
    return "Input block set";
}

}

