/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLAST_RUNNER_HPP_
#define BR_BLAST_RUNNER_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"
#include "FileWriter.hpp"

namespace npge {

/** Run blast all-against-all for given file.
\warning Output file name must be set using set_output_file() or set_rand_name()
*/
class BlastRunner : public Processor {
public:
    /** Constructor */
    BlastRunner();

protected:
    void run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
    FileWriter file_writer_;
};

}

#endif

