/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_META_PIPE_HPP_
#define NPGE_META_PIPE_HPP_

#include <string>

#include "global.hpp"

namespace npge {

/** Pipe, produced from script.
\param script Text description of pipe
\param meta Object of class Meta, used to get processors.
    If meta == 0, thread-specific default meta is used.
\param tail String where non-parsed tail of the string is written (if not 0).

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
 - max_loops number;
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
Pipe* create_pipe(const std::string& script,
                  const Meta* meta = 0,
                  std::string* tail = 0);

/** Pipe, produced from script.
Version with C-strings.
Tail, if not 0, is changed to point to beginning of tail.
*/
Pipe* create_pipe_c(const char* script,
                    const Meta* meta = 0,
                    const char** tail = 0);

/** Read script and return processors to be called.
Input is a sequence of pipe difinitions,
followed by instructions "run ProcessorName [--options];".
These processors are meant to be run by the program and they are returned.
Order of running procesors is preserved.
Pipes defined are added to Meta.
If several "run" or "add" commands define different default values
of same option, the firt one wins (see Processor::add_options()).

Example:
\code
# comment
pipe PipeName {
    name "Human readable name; Semicolon is allowed";
    bs blockset-name "blockset description";
    bs another-blockset-name "blockset description";
    max_loops 1;
    workers 2;
    no_options false;
    timing true;
    add AddBlocks;
    add Rest target=rest other=target;
    add Output target=rest;
};

#comment
pipe Pipe2 {
    max_loops 2; #comment
    add PipeName;
};

run Pipe2;
run Pipe2 --dump-seq=1;
\encode
*/
std::vector<Processor*> parse_script_to_processors(const std::string& script,
        Meta* meta);

/** Read script and return processors to be called.
This is the same to parse_script_to_processors, but it returns one processor.
If there is one "run" command, then this processor is returned.
If there is multiple "run" commands, then processors are added to pipe
named "Main pipe" and this pipe is returned.
If no "run ..." command was found, returns 0.
*/
Processor* parse_script(const std::string& script, Meta* meta);

}

#endif

