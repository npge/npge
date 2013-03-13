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

RowStorage::RowStorage(bool keep_alignment, const std::string& storage):
    keep_alignment_(keep_alignment), row_type_(storage)
{ }

AlignmentRow* RowStorage::create_row() const {
    BOOST_ASSERT(keep_alignment());
    if (row_type_ == "compact") {
        return new CompactAlignmentRow;
    } else {
        return new MapAlignmentRow;
    }
}

void RowStorage::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("import-alignment", po::value<bool>()->default_value(keep_alignment()),
     "import alignment (not only start and stop positions)")
    ("row-type", po::value<std::string>()->default_value(row_type()),
     "way of storing alignments in memory ('map' or 'compact')");
   ;
}

void RowStorage::apply_options_impl(const po::variables_map& vm) {
    if (vm.count("import-alignment")) {
        set_keep_alignment(vm["import-alignment"].as<bool>());
    }
    std::string row_type = vm["row-type"].as<std::string>();
    if (row_type != "map" && row_type != "compact") {
        throw Exception("'row-type' must be 'map' or 'compact'");
    }
    set_row_type(row_type);
}

}

