/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "Write.hpp"
#include "OriByMajority.hpp"
#include "Rest.hpp"
#include "UniqueNames.hpp"
#include "Output.hpp"

namespace npge {

Write::Write(const std::string& prefix) {
    add(new OriByMajority);
    add(new Rest, "target=target other=target");
    add(new UniqueNames);
    add(new Output(prefix));
    declare_bs("target", "Target blockset");
}

const char* Write::name_impl() const {
    return "Grace output";
}

}

