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
#include "read_file.hpp"
#include "process.hpp"
#include "string_arguments.hpp"
#include "to_s.hpp"
#include "reentrant_getenv.hpp"
#include "Meta.hpp"

namespace npge {

bool read_env(Meta* meta, const std::string& name) {
    AnyAs value = meta->get_opt(name);
    ASSERT_FALSE(value.empty());
    std::string env_value = reentrant_getenv(name);
    if (env_value.empty()) {
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

static void read_config(Meta* meta, const std::string& fname) {
    if (fname.empty()) {
        return;
    }
    std::string script = read_file(fname);
    if (script.empty()) {
        return;
    }
    StringToArgv args;
    bool debug = false;
    execute_script(script, ":cerr", args.argc(),
                   args.argv(), meta, debug);
}

void read_config(Meta* m) {
    for (int i = 0; i <= 9; i++) {
        std::string c0 = "CONFIG" + TO_S(i);
        std::string c = m->get_opt(c0, std::string()).to_s();
        if (c == "ENV") {
            read_all_env(m);
        } else if (c == "LOCAL_CONF") {
            std::string l = m->get_opt("LOCAL_CONF",
                                       std::string()).to_s();
            read_config(m, l);
        } else {
            read_config(m, c);
        }
    }
}

}

