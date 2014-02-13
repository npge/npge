/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_STATS_HPP_
#define BR_STATS_HPP_

#include "Processor.hpp"
#include "FileWriter.hpp"

namespace bloomrepeats {

/** Print human readable summary and statistics about block set */
class Stats : public Processor {
public:
    /** Constructor */
    Stats();

    /** Access file writer */
    FileWriter& file_writer() {
        return file_writer_;
    }

protected:
    bool run_impl() const;

    const char* name_impl() const;

private:
    FileWriter file_writer_;
};

}

#endif

