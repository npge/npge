/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include "boost-xtime.hpp"
#include <boost/thread/tss.hpp>
#include <boost/thread/mutex.hpp>

#include "luabind-error.hpp"
#include "luabind-format-signature.hpp"
#include <luabind/luabind.hpp>

#include "Meta.hpp"
#include "throw_assert.hpp"
#include "Exception.hpp"
#include "Processor.hpp"
#include "meta_lib.hpp"
#include "opts_lib.hpp"
#include "util_lua.hpp"
#include "model_lua.hpp"
#include "algo_lua.hpp"
#include "lua_lib.hpp"
#include "read_config.hpp"
#include "global.hpp"

namespace npge {

static void do_nothing(Meta*) {
}

static boost::thread_specific_ptr<Meta> tss_meta_(do_nothing);

MetaThreadKeeper::MetaThreadKeeper(Meta* meta) {
    prev_ = tss_meta_.get();
    tss_meta_.reset(meta);
}

MetaThreadKeeper::~MetaThreadKeeper() {
    tss_meta_.reset(prev_);
}

typedef Processor* ProcessorPtr;
typedef boost::function<ProcessorPtr()> ProcessorReturner;
typedef std::map<std::string, ProcessorReturner> ReturnerMap;

struct GlobalOption {
    Meta::AnyReturner f;
    std::string description;
    std::string section;
};

typedef std::map<std::string, GlobalOption> AnyMap;

struct MetaImpl {
    boost::thread_specific_ptr<lua_State> l_;
    // L is initialized before other members
    // L is deleted after other members
    ReturnerMap map_;
    AnyMap opts_;
    Processor* placeholder_processor_;
    MetaThreadKeeper keeper_;
    boost::mutex l_mutex_;
    bool first_l_;

    MetaImpl(Meta* meta):
        l_(&lua_close),
        placeholder_processor_(0),
        keeper_(meta),
        first_l_(true) {
    }

    ~MetaImpl() {
    }
};

struct Meta::Impl : public MetaImpl {
    Impl(Meta* meta):
        MetaImpl(meta) {
    }
};

Meta::Meta():
    impl_(new Impl(this)) {
    reset_placeholder_processor();
    add_opts(this);
    add_meta_lib(this);
    // create Lua state
    L();
    impl_->first_l_ = false;
}

Meta::~Meta() {
    delete impl_->placeholder_processor_;
    delete impl_;
}

bool Meta::has(const std::string& key) const {
    return impl_->map_.find(key) != impl_->map_.end();
}

Processor* Meta::get_plain(const std::string& key) const {
    ReturnerMap::const_iterator it = impl_->map_.find(key);
    if (it == impl_->map_.end()) {
        throw Exception("No such processor: " + key);
    }
    const ProcessorReturner& returner = it->second;
    Processor* processor = returner();
    ASSERT_TRUE(processor);
    ASSERT_EQ(processor->key(), key);
    processor->set_meta(const_cast<Meta*>(this));
    return processor;
}

SharedProcessor Meta::get(const std::string& key) const {
    return SharedProcessor(get_plain(key));
}

static std::string get_key_and_delete(const Processor* p) {
    std::string key = p->key();
    delete p;
    return key;
}

void Meta::set_returner(const ProcessorReturner& function,
                        std::string key, bool overwrite) {
    if (key.empty()) {
        Processor* p = function();
        key = get_key_and_delete(p);
    }
    ReturnerMap::iterator it = impl_->map_.find(key);
    if (it != impl_->map_.end()) {
        if (overwrite) {
            it->second = function;
        }
    } else {
        impl_->map_[key] = function;
    }
}

Strings Meta::keys() const {
    Strings result;
    BOOST_FOREACH (const ReturnerMap::value_type& kv,
                  impl_->map_) {
        result.push_back(kv.first);
    }
    return result;
}

bool Meta::empty() const {
    return impl_->map_.empty();
}

void Meta::clear() {
    impl_->map_.clear();
    impl_->opts_.clear();
}

Processor* Meta::placeholder_processor() const {
    return impl_->placeholder_processor_;
}

void Meta::reset_placeholder_processor() {
    delete impl_->placeholder_processor_;
    impl_->placeholder_processor_ = new Processor;
    impl_->placeholder_processor_->set_meta(this);
}

AnyAs Meta::get_opt(const std::string& key, const AnyAs& dflt) const {
    AnyMap::const_iterator it = impl_->opts_.find(key);
    if (it == impl_->opts_.end()) {
        return dflt;
    } else {
        const AnyReturner& f = it->second.f;
        return f();
    }
}

static AnyAs any_returner(AnyAs value) {
    return value;
}

const std::string& Meta::get_description(const std::string& k,
        const std::string& dflt) const {
    AnyMap::const_iterator it = impl_->opts_.find(k);
    if (it == impl_->opts_.end()) {
        return dflt;
    } else {
        return it->second.description;
    }
}

void Meta::set_description(const std::string& key,
                           const std::string& description) {
    impl_->opts_[key].description = description;
}

const std::string& Meta::get_section(
    const std::string& key) const {
    AnyMap::const_iterator it = impl_->opts_.find(key);
    ASSERT_TRUE(it != impl_->opts_.end());
    return it->second.section;
}

void Meta::set_section(const std::string& key,
                       const std::string& section) {
    impl_->opts_[key].section = section;
}

void Meta::set_opt(const std::string& key, const AnyAs& value,
                   const std::string& description) {
    set_opt_func(key, boost::bind(any_returner, value));
    if (!description.empty()) {
        set_description(key, description);
    }
}

void Meta::set_opt_func(const std::string& key,
                        const AnyReturner& f) {
    impl_->opts_[key].f = f;
}

bool Meta::has_opt(const std::string& key) const {
    return impl_->opts_.find(key) != impl_->opts_.end();
}

Strings Meta::opts() const {
    Strings result;
    BOOST_FOREACH (const AnyMap::value_type& kv, impl_->opts_) {
        result.push_back(kv.first);
    }
    return result;
}

Strings Meta::sections() const {
    std::set<std::string> result;
    BOOST_FOREACH (const AnyMap::value_type& k_and_v,
                  impl_->opts_) {
        result.insert(k_and_v.second.section);
    }
    return Strings(result.begin(), result.end());
}

Strings Meta::opts_of_section(
    const std::string& section) const {
    Strings result;
    BOOST_FOREACH (const AnyMap::value_type& k_and_v,
                  impl_->opts_) {
        if (k_and_v.second.section == section) {
            result.push_back(k_and_v.first);
        }
    }
    return result;
}

void Meta::remove_opt(const std::string& key) {
    impl_->opts_.erase(key);
}

void Meta::attach_to_lua(lua_State* L) {
    boost::mutex::scoped_lock lock(impl_->l_mutex_);
    impl_->l_.reset(L);
    init_util_lua(L);
    init_model_lua(L);
    init_algo_lua(L);
    using namespace luabind;
    globals(L)["meta"] = this;
    luaL_openlibs(L);
    add_lua_lib(this);
    {
        AnyMap opts = impl_->opts_;
        read_config(this);
        if (!impl_->first_l_) {
            // Restore options state.
            // Options can be changed in script or
            // in terminal.
            // This changes should be overwritten
            // while re-reading configs from
            // a worker thread.
            impl_->opts_ = opts;
        }
    }
}

lua_State* Meta::L() {
    lua_State* L = impl_->l_.get();
    if (!L) {
        L = luaL_newstate();
        attach_to_lua(L);
    }
    return L;
}

Meta* Meta::instance() {
    return tss_meta_.get();
}

}

