/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/bind.hpp>

#include "Meta.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "Processor.hpp"
#include "meta_lib.hpp"
#include "pipe_lib.hpp"
#include "opts_lib.hpp"
#include "read_config.hpp"
#include "global.hpp"

namespace npge {

Meta::Meta() {
    placeholder_processor_ = new Processor;
    placeholder_processor_->set_meta(this);
    add_opts(this);
    read_all_env(this);
    add_meta_lib(this);
    add_pipe_lib(this);
}

Meta::~Meta() {
    delete placeholder_processor_;
}

bool Meta::has(const std::string& key) const {
    return map_.find(key) != map_.end();
}

Processor* Meta::get_plain(const std::string& key) const {
    ReturnerMap::const_iterator it = map_.find(key);
    if (it == map_.end()) {
        throw Exception("No such proessor: " + key);
    }
    const ProcessorReturner& returner = it->second;
    Processor* processor = returner();
    ASSERT_EQ(processor->key(), key);
    processor->set_meta(const_cast<Meta*>(this));
    return processor;
}

SharedProcessor Meta::get(const std::string& key) const {
    return SharedProcessor(get_plain(key));
}

Strings Meta::keys() const {
    Strings result;
    BOOST_FOREACH (const ReturnerMap::value_type& key_and_func, map_) {
        result.push_back(key_and_func.first);
    }
    return result;
}

bool Meta::empty() const {
    return map_.empty();
}

void Meta::clear() {
    map_.clear();
}

AnyAs Meta::get_opt(const std::string& key, const AnyAs& dflt) const {
    AnyMap::const_iterator it = opts_.find(key);
    if (it == opts_.end()) {
        return dflt;
    } else {
        return it->second();
    }
}

static AnyAs any_returner(AnyAs value) {
    return value;
}

void Meta::set_opt(const std::string& key, const AnyAs& value) {
    set_opt_func(key, boost::bind(any_returner, value));
}

void Meta::set_opt_func(const std::string& key,
                        const AnyReturner& f) {
    opts_[key] = f;
}

Strings Meta::opts() const {
    Strings result;
    BOOST_FOREACH (const AnyMap::value_type& key_and_value, opts_) {
        result.push_back(key_and_value.first);
    }
    return result;
}

void Meta::remove_opt(const std::string& key) {
    opts_.erase(key);
}

std::string Meta::get_key_and_delete(const Processor* p) {
    std::string key = p->key();
    delete p;
    return key;
}

}

