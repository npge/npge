/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PROCESS_HPP_
#define BR_PROCESS_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Wrap a processor as a program.
\param argc Number of arguments
\param argv Arguments
\param processor Processor to be applied
\param name Name of the program
\param positional Name of positional option

Add and apply program options to the processor,
apply it to an empty block set.

Return non-zerro on error.
*/
int process(int argc, char** argv,
            Processor* processor,
            const std::string& name = "",
            const std::string& positional = "");

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

