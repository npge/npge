/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_EXTERNAL_ALIGNER_HPP_
#define BR_EXTERNAL_ALIGNER_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Align blocks with external alignment tool */
class ExternalAligner : public Processor {
public:
    /** Constructor
    \param cmd Command template. Use %1% as input of aligner, %2% as output.
    */
    ExternalAligner(const std::string& cmd);

    /** Get command template */
    const std::string& cmd() const {
        return cmd_;
    }

    /** Set command template */
    void set_cmd(const std::string& cmd) {
        cmd_ = cmd;
    }

    /** Apply external aligner to a blick */
    void align_block(Block* block) const;

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;

private:
    std::string cmd_;
};

}

#endif

