/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <typeinfo>
#include <iostream>
#include <boost/thread.hpp>
#include <boost/thread/tss.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#include "Processor.hpp"
#include "BlockSet.hpp"

namespace bloomrepeats {

Processor::Processor():
    workers_(1), no_options_(false), timing_(false), milliseconds_(0)
{ }

Processor::~Processor() {
    if (timing()) {
        using namespace boost::posix_time;
        std::cerr << name_ << ": ";
        std::cerr << to_simple_string(milliseconds(milliseconds_)) << std::endl;
    }
}

BlockSetPtr Processor::block_set() const {
    return target_block_set_.other();
}

void Processor::set_block_set(BlockSetPtr block_set) {
    target_block_set_.set_other(block_set);
}

void Processor::set_target_processor(Processor* processor) {
    target_block_set_.set_processor(processor);
}

void Processor::set_target_other(OtherBlockSet* other) {
    target_block_set_.set_other_block_set(other);
}

void Processor::set_empty_block_set() {
    set_block_set(boost::make_shared<BlockSet>());
}

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
         "number of threads")
        ("timing", "measure time for each processor")
       ;
        if (recursive_options()) {
            add_options_impl(desc);
        } else {
            po::options_description temp;
            add_options_impl(temp);
            po::options_description new_opts(name());
            add_new_options(temp, new_opts, &desc);
            if (!new_opts.options().empty()) {
                desc.add(new_opts);
            }
        }
    }
}

void Processor::apply_options(const po::variables_map& vm) {
    if (!no_options()) {
        apply_options_impl(vm);
        set_workers(vm["workers"].as<int>());
    }
    if (vm.count("timing")) {
        set_timing(true);
    }
}

bool Processor::run() const {
    bool result = false;
    if (workers() != 0 && block_set()) {
        using namespace boost::posix_time;
        ptime before, after;
        if (timing()) {
            before = microsec_clock::universal_time();
        }
        result = run_impl();
        if (timing()) {
            after = microsec_clock::universal_time();
            milliseconds_ += (after - before).total_milliseconds();
            name_ = name();
            if (name_.empty()) {
                name_ = typeid(*this).name();
            }
        }
    }
    return result;
}

const char* Processor::name() const {
    return name_impl();
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

const char* Processor::name_impl() const {
    return "";
}

static boost::thread_specific_ptr<bool> flag_;

static bool flag() {
    if (flag_.get() == 0) {
        flag_.reset(new bool(false));
    }
    return *flag_;
}

static void set_flag(bool value) {
    if (flag_.get() == 0) {
        flag_.reset(new bool(false));
    }
    *flag_ = value;
}

bool Processor::recursive_options() const {
    po::options_description temp;
    set_flag(false);
    add_options_impl(temp);
    bool result = flag() == true;
    set_flag(true);
    return result;
}

}

