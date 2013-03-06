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
#include "SizeLimits.hpp"

namespace bloomrepeats {

/** Print human readable summary and statistics about block set */
class Stats : public Processor, public FileWriter, public SizeLimits {
protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;
};

}

#endif

