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
    OptionsPrefix(prefix),
    processor_(processor),
    opts_(opts),
    p_(0)
{ }

void MetaProcessor::add_options_impl(po::options_description& desc) const {
    add_unique_options(desc)
    ("processor", po::value<std::string>()->default_value(processor()),
     "processor (implementation) name")
    ("opts", po::value<std::string>()->default_value(opts()),
     "options for processor (implementation)")
   ;
}

void MetaProcessor::apply_options_impl(const po::variables_map& vm) {
    if (vm.count(prefixed("processor"))) {
        set_processor(vm[prefixed("processor")].as<std::string>());
    }
    if (vm.count(prefixed("opts"))) {
        set_opts(vm[prefixed("opts")].as<std::string>());
    }
}

bool MetaProcessor::run_impl() const {
    if (!meta()->has(processor())) {
        throw Exception("No processor '" + processor() +
                        "' found for MetaProcessor.");
    }
    if (!p_) {
        p_ = meta()->get_plain(processor());
        p_->set_parent(const_cast<MetaProcessor*>(this)); // FIXME
        p_->set_workers(workers());
        p_->set_timing(timing());
        p_->set_options(opts(), const_cast<MetaProcessor*>(this)); // FIXME
    }
    return p_->run();
}

const char* MetaProcessor::name_impl() const {
    return "Any processor";
}

}

