/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "pipe_lib.hpp"
#include "meta_pipe.hpp"
#include "Meta.hpp"
#include "Pipe.hpp"
#include "Processor.hpp"

namespace npge {

void add_pipe_c(Meta* meta, const char* script,
                int length) {
    const char** z = 0;
    meta->set_returner(boost::bind(create_pipe_c, script,
                                   meta, z, length));
}

void add_pipe(Meta* meta, const std::string& script) {
    std::string* z = 0;
    meta->set_returner(boost::bind(create_pipe, script, meta, z));
}

#define NPGE_SCRIPT(...) #__VA_ARGS__

void add_pipe_lib(Meta* meta) {
    add_pipe_c(meta, NPGE_SCRIPT(
    pipe Test {
        add Hash;
    };
             ));
}

}

