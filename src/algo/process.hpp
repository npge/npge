/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
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

/** Wrap a processor as a program.
\param argc Number of arguments
\param argv Arguments
\param processor Processor to be applied
\param name Name of the program
\param positional Name of positional option
\param catch_sigint If SIGINT should be catched by the program

Add and apply program options to the processor,
apply it to an empty block set.

Return non-zerro on error.
*/
int process(int argc, char** argv,
            Processor* processor,
            const std::string& name = "",
            const std::string& positional = "",
            bool catch_sigint = true,
            bool print_changed = false);

/** Process each processor (until first error) and delete them all */
int process_and_delete(int argc, char** argv,
                       const std::vector<Processor*>& processors,
                       const std::string& positional = "",
                       bool print_changed = false);

/** Run commands from script */
int execute_script(const std::string& script,
                   const std::string& output,
                   int argc, char** argv, Meta* meta,
                   bool debug = false,
                   const std::string& positional = "in-blocks",
                   bool print_changed = false);

/** Run commands from input file one by one.
See name_to_istream().
Commands should be in format parse_script_to_processors() accepts.
Return 0 if no error occured.
Return code of last error is any. For parsing error code 15 is used.
*/
int interactive_loop(const std::string& input, const std::string& output,
                     int argc, char** argv, Meta* meta);

/** Applies options of source to destination.
Unknown options are ignored.
*/
void copy_processor_options(Processor& dest, const Processor& source);

}

#endif

