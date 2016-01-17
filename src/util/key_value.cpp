/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "key_value.hpp"

namespace npge {

std::string extract_value(const std::string& values, const std::string& key) {
    using boost::tokenizer;
    using boost::escaped_list_separator;
    typedef tokenizer<escaped_list_separator<char> > tok_t;
    tok_t t(values, escaped_list_separator<char>('\\', ' ', '\"'));
    std::string look_for = key + "=";
    BOOST_FOREACH (std::string opt, t) {
        using namespace boost::algorithm;
        if (starts_with(opt, look_for)) {
            if (opt.length() > look_for.length()) {
                return opt.substr(look_for.length());
            } else {
                return "";
            }
        }
    }
    return "";
}

}

