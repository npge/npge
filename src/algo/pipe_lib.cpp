/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <memory>
#include <boost/foreach.hpp>

#include "pipe_lib.hpp"
#include "meta_pipe.hpp"
#include "Meta.hpp"
#include "Pipe.hpp"
#include "Processor.hpp"

namespace npge {

void add_pipe_c(Meta* meta,
                const char* script,
                const std::string& key,
                int length) {
    const char** z = 0;
    meta->set_returner(boost::bind(create_pipe_c, script,
                                   meta, z, length), key);
}

void add_pipe(Meta* meta, const std::string& script,
              const std::string& key) {
    std::string* z = 0;
    meta->set_returner(boost::bind(create_pipe,
                                   script, meta, z), key);
}

void add_pipes_c(Meta* m, const char* s) {
    int length = strlen(s);
    const char* tail;
    while (true) {
        std::auto_ptr<Pipe> p;
        try {
            p.reset(create_pipe_c(s, m, &tail, length));
        } catch (...) {
            break;
        }
        int l = tail - s;
        add_pipe_c(m, s, p->key(), l);
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

void add_pipe_lib(Meta* meta) {
#include "pipe_lib.npge"
}

}

