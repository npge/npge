/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/bind.hpp>

#include "RowStorage.hpp"
#include "Processor.hpp"
#include "AlignmentRow.hpp"

namespace npge {

static bool check_row_type(std::string& message, Processor* p) {
    if (p->opt_value("row-type").as<std::string>() != "map" &&
            p->opt_value("row-type").as<std::string>() != "compact") {
        message = "row-type must be 'map' or 'compact'";
        return false;
    }
    return true;
}

void add_row_storage_options(Processor* p) {
    p->add_opt("row-type",
               "way of storing alignments in memory ('map' or 'compact')",
               std::string("compact"));
    p->add_opt_check(boost::bind(check_row_type, _1, p));
}

RowType row_type(const Processor* p) {
    return (p->opt_value("row-type").as<std::string>() == "map") ?
           MAP_ROW : COMPACT_ROW;
}

AlignmentRow* create_row(const Processor* p) {
    return AlignmentRow::new_row(row_type(p));
}

}

