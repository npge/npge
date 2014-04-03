/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PIPE_LIB_HPP_
#define BR_PIPE_LIB_HPP_

#include "global.hpp"

namespace bloomrepeats {

/** Add pipe as script.
Const C-string must remain available until deletion of Meta.
*/
void add_pipe(Meta* meta, const char* script);

/** Add pipe as script */
void add_pipe(Meta* meta, const std::string& script);

/** Add standard pipes collection */
void add_pipe_lib(Meta* meta);

}

#endif

