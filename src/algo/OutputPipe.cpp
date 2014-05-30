/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "OutputPipe.hpp"
#include "OriByMajority.hpp"
#include "Connector.hpp"
#include "Rest.hpp"
#include "UniqueNames.hpp"
#include "Output.hpp"

namespace npge {

OutputPipe::OutputPipe(const std::string& prefix) {
    add(new OriByMajority);
    add(new Connector);
    add(new Rest, "target=target other=target");
    add(new Connector);
    add(new UniqueNames);
    add(new Output(prefix));
    declare_bs("target", "Target blockset");
}

const char* OutputPipe::name_impl() const {
    return "Grace output";
}

}

