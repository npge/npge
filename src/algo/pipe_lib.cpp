/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "pipe_lib.hpp"
#include "meta_pipe.hpp"
#include "Meta.hpp"
#include "Pipe.hpp"
#include "Processor.hpp"

namespace bloomrepeats {

void add_pipe(Meta* meta, const char* script) {
    const char** z = 0;
    meta->set_returner(boost::bind(create_pipe_c, script, meta, z));
}

void add_pipe(Meta* meta, const std::string& script) {
    std::string* z = 0;
    meta->set_returner(boost::bind(create_pipe, script, meta, z));
}

#define BR_SCRIPT(...) #__VA_ARGS__

void add_pipe_lib(Meta* meta) {
    add_pipe(meta, BR_SCRIPT(
    pipe Test {
        add Hash;
    };
             ));
}

}
