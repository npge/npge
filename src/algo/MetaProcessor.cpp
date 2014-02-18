/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "MetaProcessor.hpp"
#include "Meta.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

MetaProcessor::MetaProcessor(const std::string& prefix,
                             const std::string& processor,
                             const std::string& opts):
    p_(0) {
    set_opt_prefix(prefix);
    add_opt("processor", "processor name", processor);
    add_opt("opts", "options for processor", opts);
}

bool MetaProcessor::run_impl() const {
    std::string processor = opt_value("processor").as<std::string>();
    std::string opts = opt_value("opts").as<std::string>();
    if (!meta()->has(processor)) {
        throw Exception("No processor '" + processor +
                        "' found for MetaProcessor.");
    }
    if (!p_) {
        p_ = meta()->get_plain(processor);
        p_->set_parent(const_cast<MetaProcessor*>(this)); // FIXME
        p_->set_workers(workers());
        p_->set_timing(timing());
        p_->set_options(opts, const_cast<MetaProcessor*>(this)); // FIXME
    }
    return p_->run();
}

const char* MetaProcessor::name_impl() const {
    return "Any processor";
}

}

