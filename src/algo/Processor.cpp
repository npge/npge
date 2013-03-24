/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <typeinfo>
#include <boost/thread.hpp>
#include <boost/thread/tss.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#include "Processor.hpp"
#include "class_name.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"

namespace bloomrepeats {

struct BlockSetHolder {
public:
    BlockSetHolder():
        processor_(0)
    { }

    BlockSetPtr block_set() const {
        if (block_set_) {
            return block_set_;
        } else if (processor_) {
            return processor_->get_bs(name_);
        } else {
            block_set_ = new_bs();
            return block_set_;
        }
    }

    void set_block_set(BlockSetPtr block_set) {
        block_set_ = block_set;
        processor_ = 0;
        name_.clear();
    }

    void set_processor(Processor* processor, const std::string& name) {
        block_set_.reset();
        processor_ = processor;
        name_ = name;
    }

private:
    mutable BlockSetPtr block_set_;
    // or
    Processor* processor_;
    std::string name_;
};

typedef std::map<std::string, BlockSetHolder> BlockSetMap;

struct Processor::Impl {
    Impl():
        workers_(1), no_options_(false), timing_(false), milliseconds_(0)
    { }

    BlockSetMap map_;
    int workers_;
    bool no_options_;
    bool timing_;
    mutable int milliseconds_;
    std::string name_;
    po::options_description ignored_options_;
    std::string key_;
};

Processor::Processor() {
    impl_ = new Impl;
}

Processor::~Processor() {
    if (timing()) {
        using namespace boost::posix_time;
        std::cerr << key() << ": ";
        std::cerr << to_simple_string(milliseconds(impl_->milliseconds_));
        std::cerr << std::endl;
    }
    delete impl_;
}

BlockSetPtr Processor::get_bs(const std::string& name) const {
    BlockSetMap::const_iterator it = impl_->map_.find(name);
    if (it == impl_->map_.end()) {
        BlockSetPtr bs = new_bs();
        impl_->map_[name].set_block_set(bs);
        return bs;
    } else {
        return it->second.block_set();
    }
}

void Processor::set_bs(const std::string& name, BlockSetPtr bs) {
    impl_->map_[name].set_block_set(bs);
}

void Processor::point_bs(const std::string& mapping, Processor* processor) {
    size_t eq_pos = mapping.find("=");
    BOOST_ASSERT_MSG(eq_pos != std::string::npos,
                     ("Bad mapping: " + mapping).c_str());
    std::string name_in_this = mapping.substr(0, eq_pos);
    std::string name_in_processor = mapping.substr(eq_pos + 1);
    impl_->map_[name_in_this].set_processor(processor, name_in_processor);
}

BlockSetPtr Processor::block_set() const {
    return get_bs("target");
}

void Processor::set_block_set(BlockSetPtr block_set) {
    set_bs("target", block_set);
}

BlockSetPtr Processor::other() const {
    return get_bs("other");
}

void Processor::set_other(BlockSetPtr other) {
    set_bs("other", other);
}

void Processor::set_empty_block_set() {
    set_block_set(new_bs());
}

void Processor::set_empty_other() {
    set_other(new_bs());
}

int Processor::workers() const {
    return impl_->workers_;
}

void Processor::set_workers(int workers) {
    if (workers == -1) {
        impl_->workers_ = boost::thread::hardware_concurrency();
        if (workers == 0) {
            impl_->workers_ = 1;
        }
    }
    impl_->workers_ = workers;
}

bool Processor::no_options() const {
    return impl_->no_options_;
}

void Processor::set_no_options(bool no_options) {
    impl_->no_options_ = no_options;
}

void Processor::add_ignored_option(const std::string& option) {
    add_unique_options(impl_->ignored_options_)(option.c_str(), "");
}

bool Processor::timing() const {
    return impl_->timing_;
}

void Processor::set_timing(bool timing) {
    impl_->timing_ = timing;
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
            po::options_description not_ignored;
            add_new_options(temp, not_ignored, &impl_->ignored_options_);
            po::options_description new_opts(name());
            add_new_options(not_ignored, new_opts, &desc);
            if (!new_opts.options().empty()) {
                desc.add(new_opts);
            }
        }
    }
}

void Processor::apply_options(const po::variables_map& vm) {
    if (!no_options()) {
        apply_options_impl(vm);
        if (vm.count("workers")) {
            set_workers(vm["workers"].as<int>());
            if (std::abs(vm["workers"].as<int>()) < 1) {
                throw Exception("'workers' number must be >= 1");
            }
        }
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
            impl_->milliseconds_ += (after - before).total_milliseconds();
        }
    }
    return result;
}

std::string Processor::name() const {
    if (!impl_->name_.empty()) {
        return impl_->name_.c_str();
    } else {
        const char* ni = name_impl();
        if (*ni == '\0') {
            // empty
            return key();
        } else {
            return ni;
        }
    }
}

void Processor::set_name(const std::string& name) {
    impl_->name_ = name;
}

bool Processor::apply(const BlockSetPtr& bs) const {
    BlockSetPtr prev = block_set();
    const_cast<Processor*>(this)->set_block_set(bs);
    bool result = run();
    const_cast<Processor*>(this)->set_block_set(prev);
    return result;
}

std::string Processor::key() const {
    if (impl_->key_.empty()) {
        return processor_name(this);
    } else {
        return impl_->key_;
    }
}

void Processor::set_key(const std::string& key) {
    impl_->key_ = key;
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

std::string processor_name(const Processor* processor) {
    return class_name(typeid(*processor).name());
}

std::string processor_name(const ProcessorPtr& processor) {
    return processor_name(processor.get());
}

}

