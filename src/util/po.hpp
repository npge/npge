/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PO_HPP_
#define BR_PO_HPP_

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

#include "global.hpp"

namespace bloomrepeats {

#ifndef DOXYGEN_ONLY
class AddUniqueOptions : private po::options_description,
    public po::options_description_easy_init {
public:
    AddUniqueOptions(po::options_description& desc);
    ~AddUniqueOptions();
private:
    po::options_description& desc_;
};
#endif

/** Add only new options.
Usage:
\code
add_unique_options(desc)
    ("help,h", "produce help message")
;
\endcode
*/
AddUniqueOptions add_unique_options(po::options_description& desc);

/** Add general options like "--help" */
void add_general_options(po::options_description& desc);

/** Read options. Return non-zerro on error */
int read_options(int argc, char** argv, po::variables_map& vm,
                 const po::options_description& desc,
                 const po::positional_options_description& pod);

}

#endif

