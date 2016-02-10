/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "MkDir.hpp"
#include "name_to_stream.hpp"

namespace npge {

MkDir::MkDir() {
    add_opt("dirname", "Dir to be created",
            std::string(), true);
}

void MkDir::run_impl() const {
    make_dir(opt_value("dirname").as<std::string>());
}

const char* MkDir::name_impl() const {
    return "Create directory";
}

}

