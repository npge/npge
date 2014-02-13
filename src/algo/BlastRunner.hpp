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
class BlastRunner : public Processor {
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

    /** Return if low complexity regions are skipped */
    bool skip_low_complexity_regions() const {
        return skip_low_complexity_regions_;
    }

    /** Set if low complexity regions are skipped.
    Defaults to false (-F F).
    */
    void set_skip_low_complexity_regions(bool value) {
        skip_low_complexity_regions_ = value;
    }

protected:
    bool run_impl() const;

    const char* name_impl() const;

private:
    FileReader file_reader_;
    FileWriter file_writer_;
    float evalue_;
    bool skip_low_complexity_regions_;
};

}

#endif

