/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_META_PIPE_HPP_
#define BR_META_PIPE_HPP_

#include <string>

#include "global.hpp"

namespace bloomrepeats {

/** Pipe, produced from script.
\param script Text description of pipe
\param meta Object of class Meta, used to get processors.
    If meta == 0, default meta is used.

Comments are started from '#', until end of line.
Pipe description starts with "pipe".
Pipe key is next token. Pipe can be unnamed.
Rest of pipe description is placed in { }.
Each block is ended with ';'.
Commands:
 - name "name of this pipe";
 - workers number-of-workers;
 - no_options true/false;
 - timing true/false;
 - add ProcessorName options;

Example:
\code
# comment
pipe PipeName {
    name "Human readable name; Semicolon is allowed";
    max_loops 1;
    workers 2;
    no_options false;
    timing true;
    add AddBlocks;
    add Rest target=rest other=target;
    add Output target=rest;
};
\encode
*/
boost::shared_ptr<Pipe> create_pipe(const std::string& script,
        const Meta* meta = 0);

}

#endif

