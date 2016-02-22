/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_REPORT_LIST_HPP_
#define NPGE_REPORT_LIST_HPP_

#include <ostream>
#include <algorithm>
#include <boost/format.hpp>

#include "boundaries.hpp"

namespace npge {

/** Print properties of list of pos_t or double */
template<typename Vector>
void report_list(
    std::ostream& o,
    const Vector& list,
    const char* prefix = "  "
) {
    if (!list.empty()) {
        typename Vector::value_type min_value, max_value;
        min_value = *std::min_element(list.begin(), list.end());
        max_value = *std::max_element(list.begin(), list.end());
        if (min_value == max_value) {
            o << prefix << "all:\t" << min_value << "\n";
        } else {
            o << prefix << "min:\t" << min_value << "\n";
            o << prefix << "median:\t" << median_element(list) << "\n";
            std::streamsize precision = o.precision();
            o.setf(std::ios::fixed, std::ios_base::floatfield);
            o.precision(2);
            o << prefix << "avg:\t" << avg_element_double(list) << "\n";
            o.precision(precision);
            o.unsetf(std::ios_base::floatfield);
            o << prefix << "max:\t" << max_value << "\n";
        }
    }
    o << std::endl;
}

}

#endif

