/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
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

/** Print properties of list of size_t or float */
template<typename Vector>
static void report_list(std::ostream& o, const Vector& list) {
    o << " number=" << list.size();
    if (!list.empty()) {
        typename Vector::value_type min_value, max_value;
        min_value = *std::min_element(list.begin(), list.end());
        max_value = *std::max_element(list.begin(), list.end());
        if (min_value == max_value) {
            o << " all=" << min_value;
        } else {
            o << " min=" << min_value;
            o << " median=" << median_element(list);
            double avg = avg_element_double(list);
            boost::format double_2("%.2f");
            o << " avg=" << str(double_2 % avg);
            o << " max=" << max_value;
        }
    }
    o << std::endl;
}

}

#endif

