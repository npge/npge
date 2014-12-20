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
#include "RawWrite.hpp"

namespace npge {

Write::Write(const std::string& prefix) {
    add(new OriByMajority);
    add(new Rest, "target=target other=target");
    add(new UniqueNames);
    add(new RawWrite(prefix));
    declare_bs("target", "Target blockset");
}

const char* Write::name_impl() const {
    return "Write target blockset to file in bs format. "
           "Uncovered parts are written as blocks";
}

}

