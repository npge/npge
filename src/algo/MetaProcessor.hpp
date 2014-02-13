/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_META_PROCESSOR_HPP_
#define BR_META_PROCESSOR_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

/** Any processor.
Processor is specified through an option --[prefix]processor.
Options (Processor::set_options()) for internal processor are specified
through an option --[prefix]opts.

run() will fail if no processor was set.
*/
class MetaProcessor : public Processor {
public:
    /** Constructor */
    MetaProcessor(const std::string& prefix = "",
                  const std::string& processor = "",
                  const std::string& opts = "");

    /** Return processor name */
    const std::string& processor() const {
        return processor_;
    }

    /** Set processor name */
    void set_processor(const std::string& processor) {
        processor_ = processor;
    }

    /** Return processor options */
    const std::string& opts() const {
        return opts_;
    }

    /** Set processor options */
    void set_opts(const std::string& opts) {
        opts_ = opts;
    }

protected:
    void add_options_impl(po::options_description& desc) const;

    void apply_options_impl(const po::variables_map& vm);

    bool run_impl() const;

    const char* name_impl() const;

private:
    std::string processor_;
    std::string opts_;

    mutable Processor* p_;
};

}

#endif

