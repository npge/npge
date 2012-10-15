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
            ProcessorPtr processor,
            const std::string& name = "",
            const std::string& positional = "");

/** Wrap a processor as a program.
Overloaded method.
*/
int process(int argc, char** argv,
            Processor* processor,
            const std::string& name = "",
            const std::string& positional = "");

}

#endif

