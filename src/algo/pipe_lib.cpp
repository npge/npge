/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

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

void add_pipes_c(Meta* m, const char* s) {
    int length = strlen(s);
    const char* tail;
    while (true) {
        Pipe* p;
        try {
            p = create_pipe_c(s, m, &tail, length);
        } catch (...) {
            break;
        }
        delete p;
        int l = tail - s;
        add_pipe_c(m, s, l);
        s = tail;
        length -= l;
    }
}

void add_pipes(Meta* meta, const std::string& script) {
    typedef std::vector<Processor*> Processors;
    Processors pp = parse_script_to_processors(script, meta);
    BOOST_FOREACH (Processor* p, pp) {
        delete p;
    }
}

#define NPGE_SCRIPT(...) #__VA_ARGS__

void add_pipe_lib(Meta* meta) {
    add_pipes_c(meta, NPGE_SCRIPT(
    pipe Test {
        add Hash;
    };
             ));
}

}

