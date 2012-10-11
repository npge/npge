/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Pipe.hpp"

namespace bloomrepeats {

Pipe::Pipe():
    max_loops_(1)
{ }

Pipe& Pipe::add(const ProcessorPtr& processor) {
    processors_.push_back(processor);
    return *this;
}

Pipe& Pipe::add(Processor* processor) {
    add(ProcessorPtr(processor));
    return *this;
}

void Pipe::add_options_impl(po::options_description& desc) const {
    BOOST_FOREACH (const ProcessorPtr& processor, processors_) {
        processor->add_options(desc);
    }
}

void Pipe::apply_options_impl(const po::variables_map& vm) {
    BOOST_FOREACH (const ProcessorPtr& processor, processors_) {
        processor->apply_options(vm);
    }
}

bool Pipe::run_impl() const {
    BOOST_FOREACH (const ProcessorPtr& processor, processors_) {
        processor->set_workers(workers());
        processor->set_block_set(block_set());
    }
    bool result = false;
    for (int i = 0; i < max_loops() || max_loops() == -1; i++) {
        bool changed = false;
        BOOST_FOREACH (const ProcessorPtr& processor, processors_) {
            changed |= processor->run();
        }
        result |= changed;
        if (!changed) {
            break;
        }
    }
    return result;
}

}

