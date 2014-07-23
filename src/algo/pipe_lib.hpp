/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PIPE_LIB_HPP_
#define NPGE_PIPE_LIB_HPP_

#include "global.hpp"

namespace npge {

/** Add pipe as script.
Const C-string must remain available until deletion of Meta.
If length = 0, strlen(script) is used as length.
If key is empty, pipe is generated to get its key.
*/
void add_pipe_c(Meta* meta,
                const char* script,
                const std::string& key = "",
                int length = 0);

/** Add pipe as script.
If key is empty, pipe is generated to get its key.
*/
void add_pipe(Meta* meta, const std::string& script,
              const std::string& key = "");

/** Add multiple pipes to meta.
Const C-string must remain available until deletion of Meta.
*/
void add_pipes_c(Meta* meta, const char* script);

/** Add multiple pipes as script */
void add_pipes(Meta* meta, const std::string& script);

/** Add standard pipes collection */
void add_pipe_lib(Meta* meta);

}

#endif

