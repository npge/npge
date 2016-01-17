/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_OPTS_LIB_HPP_
#define NPGE_OPTS_LIB_HPP_

#include "global.hpp"

namespace npge {

/** Add global options to Meta */
void add_opts(Meta* meta);

/** Overwrite value of LOCAL_CONF */
void set_local_conf(const std::string& conf);

/** Make external command */
std::string make_external_cmd(const Meta* meta,
                              const std::string& name);

}

#endif

