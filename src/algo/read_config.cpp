/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdlib>
#include <boost/foreach.hpp>

#include "throw_assert.hpp"
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

}

