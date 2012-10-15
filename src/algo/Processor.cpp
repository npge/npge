/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/thread.hpp>

#include "Processor.hpp"

namespace bloomrepeats {

Processor::Processor():
    workers_(1), no_options_(false)
{ }

void Processor::set_workers(int workers) {
    if (workers == -1) {
        workers_ = boost::thread::hardware_concurrency();
        if (workers == 0) {
            workers_ = 1;
        }
    }
    workers_ = workers;
}

void Processor::assign(const Processor& other) {
    set_block_set(other.block_set());
    set_workers(other.workers());
}

void Processor::add_options(po::options_description& desc) const {
    if (!no_options()) {
        add_unique_options(desc)
        ("workers", po::value<int>()->default_value(workers()),
         "number of threads used to find anchors")
       ;
        add_options_impl(desc);
    }
}

void Processor::apply_options(const po::variables_map& vm) {
    if (!no_options()) {
        apply_options_impl(vm);
        set_workers(vm["workers"].as<int>());
    }
}

bool Processor::run() const {
    if (workers() != 0 && block_set()) {
        return run_impl();
    }
    return false;
}

bool Processor::apply(const BlockSetPtr& bs) const {
    BlockSetPtr prev = block_set();
    const_cast<Processor*>(this)->set_block_set(bs);
    bool result = run();
    const_cast<Processor*>(this)->set_block_set(prev);
    return result;
}

void Processor::add_options_impl(po::options_description& desc) const
{ }

void Processor::apply_options_impl(const po::variables_map& vm)
{ }

bool Processor::run_impl() const {
    return false;
}

}

