/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PROCESS_HPP_
#define NPGE_PROCESS_HPP_

#include "Processor.hpp"

namespace npge {

/** Print processor children as tree*/
void print_processor_tree(const std::string& output,
                          Processor* processor, int indent = 0);

/** Print help message for the processor */
void print_help(const std::string& output, const Processor* processor,
                const std::string& app = "app",
                const std::string& positional = "");

/** Print current configuration to config file */
void print_config(const std::string& out, const Meta* meta);

}

#endif

