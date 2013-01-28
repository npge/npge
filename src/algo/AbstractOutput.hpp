/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ABSTRACT_OUTPUT_HPP_
#define BR_ABSTRACT_OUTPUT_HPP_

#include <iosfwd>

#include "Processor.hpp"
#include "OptionsPrefix.hpp"

namespace bloomrepeats {

/** Print some function from blocks to file or to stdout */
class AbstractOutput : public Processor, public OptionsPrefix {
public:
    /** Get output file with all blocks */
    const std::string& file() const {
        return file_;
    }

    /** Set output file with all blocks */
    void set_file(const std::string& file) {
        file_ = file;
    }

    /** Get mask of output files (${block} is replaced with block name) */
    const std::string& mask() const {
        return mask_;
    }

    /** Set mask of output files */
    void set_mask(const std::string& mask) {
        mask_ = mask;
    }

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    /** Print block to output stream */
    virtual void print_block(std::ostream& o, Block* block) const = 0;

private:
    std::string file_;
    std::string mask_;
};

}

#endif

