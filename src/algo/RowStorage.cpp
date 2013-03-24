/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "RowStorage.hpp"
#include "AlignmentRow.hpp"
#include "Exception.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

RowStorage::RowStorage(bool keep_alignment, RowType row_type):
    keep_alignment_(keep_alignment), row_type_(row_type)
{ }

AlignmentRow* RowStorage::create_row() const {
    BOOST_ASSERT(keep_alignment());
    return AlignmentRow::new_row(row_type());
}

void RowStorage::add_options_impl(po::options_description& desc) const {
    std::string row_type_str = row_type() == COMPACT_ROW ? "compact" : "map";
    add_unique_options(desc)
    ("import-alignment", po::value<bool>()->default_value(keep_alignment()),
     "import alignment (not only start and stop positions)")
    ("row-type", po::value<std::string>()->default_value(row_type_str),
     "way of storing alignments in memory ('map' or 'compact')");
   ;
}

void RowStorage::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("import-alignment")) {
        set_keep_alignment(vm["import-alignment"].as<bool>());
    }
    if (vm.count("row-type")) {
        std::string row_type_str = vm["row-type"].as<std::string>();
        if (row_type_str == "map") {
            set_row_type(MAP_ROW);
        } else if (row_type_str == "compact") {
            set_row_type(COMPACT_ROW);
        } else {
            throw Exception("'row-type' must be 'map' or 'compact'");
        }
    }
}

}

