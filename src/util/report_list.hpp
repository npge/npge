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
void report_list(std::ostream& o, const Vector& list) {
    o << " number:\t" << list.size() << "\n";
    if (!list.empty()) {
        typename Vector::value_type min_value, max_value;
        min_value = *std::min_element(list.begin(), list.end());
        max_value = *std::max_element(list.begin(), list.end());
        if (min_value == max_value) {
            o << " all:\t" << min_value << "\n";
        } else {
            o << " min:\t" << min_value << "\n";
            o << " median:\t" << median_element(list) << "\n";
            std::streamsize precision = o.precision();
            o.setf(std::ios::fixed, std::ios_base::floatfield);
            o.precision(2);
            o << " avg:\t" << avg_element_double(list) << "\n";
            o.precision(precision);
            o.unsetf(std::ios_base::floatfield);
            o << " max:\t" << max_value << "\n";
        }
    }
    o << std::endl;
}

}

#endif

