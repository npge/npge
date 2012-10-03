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
    workers_(1)
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
    add_options_impl(desc);
    desc.add_options()
    ("workers", po::value<int>()->default_value(workers()),
     "number of threads used to find anchors")
   ;
}

void Processor::apply_options(const po::variables_map& vm) {
    apply_options_impl(vm);
    set_workers(vm["workers"].as<int>());
}

bool Processor::run() const {
    if (workers() != 0 && block_set()) {
        return run_impl();
    }
    return false;
}

void Processor::add_options_impl(po::options_description& desc) const
{ }

void Processor::apply_options_impl(const po::variables_map& vm)
{ }

bool Processor::run_impl() const {
    return false;
}

}

