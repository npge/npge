/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "throw_assert.hpp"
#include "name_to_stream.hpp"
#include "Meta.hpp"

namespace npge {

bool read_env(Meta* meta, const std::string& name) {
    AnyAs value = meta->get_opt(name);
    ASSERT_FALSE(value.empty());
    char* env_value = getenv(name.c_str());
    if (!env_value) {
        return false;
    }
    value.from_s(env_value);
    meta->set_opt(name, value);
    return true;
}

void read_all_env(Meta* meta) {
    BOOST_FOREACH (const std::string& name, meta->opts()) {
        read_env(meta, name);
    }
}

void read_config_file(Meta* meta, const std::string& cfg) {
    typedef boost::shared_ptr<std::istream> IStreamPtr;
    IStreamPtr input = name_to_istream(cfg);
    for (std::string line; std::getline(*input, line);) {
        using namespace boost::algorithm;
        Strings parts;
        split(parts, line, is_any_of("="), token_compress_on);
        ASSERT_EQ(parts.size(), 2);
        std::string& name = parts[0];
        std::string& cfg_value = parts[1];
        trim(name);
        trim(cfg_value);
        AnyAs value = meta->get_opt(name);
        if (!value.empty()) {
            value.from_s(cfg_value);
        } else {
            // read unknown options as strings
            value = cfg_value;
        }
        meta->set_opt(name, value);
    }
}

}



