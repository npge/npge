/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "Meta.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "Processor.hpp"
#include "meta_lib.hpp"
#include "pipe_lib.hpp"
#include "global.hpp"

namespace bloomrepeats {

Meta::Meta() {
    placeholder_processor_ = new Processor;
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

std::string Meta::get_key_and_delete(const Processor* p) {
    std::string key = p->key();
    delete p;
    return key;
}

}

