/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_BLAST_RUNNER_HPP_
#define BR_BLAST_RUNNER_HPP_

#include "Processor.hpp"
#include "FileReader.hpp"
#include "FileWriter.hpp"

namespace bloomrepeats {

/** Run blast all-against-all for given file.
\warning Output file name must be set using set_output_file() or set_rand_name()
*/
class BlastRunner : public Processor, public FileReader, public FileWriter {
public:
    /** Constructor */
    BlastRunner();

    /** Get e-value threshold */
    float evalue() const {
        return evalue_;
    }

    /** Set e-value threshold.
    Defaults to 0.001.
    */
    void set_evalue(float evalue) {
        evalue_ = evalue;
    }

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;

private:
    float evalue_;
};

}

#endif

