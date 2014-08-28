/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_META_HPP_
#define NPGE_META_HPP_

#include <string>
#include <map>
#include <vector>
#include <boost/utility.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include <lua.hpp>

#include "global.hpp"
#include "AnyAs.hpp"

namespace npge {

/** Return processor by key */
class Meta : boost::noncopyable {
public:
    /** Constructor.
    Add returners of all processors from algo/.
    */
    Meta();

    /** Destructor */
    virtual ~Meta();

    /** Return if a processor is associated with the key */
    bool has(const std::string& key) const;

    /** Return plain pointer to processor instance by key.
    If no processor were associated with the key,
    assertation exception is thrown.

    Processor.meta() of returned processor is set to this.
    */
    Processor* get_plain(const std::string& key) const;

    /** Return shared ptr to processor instance by key.
    \see get_plain()
    */
    SharedProcessor get(const std::string& key) const;

    /** Associate function returning processor.
    If key is empty, it is taken as function()->key().
    */
    template<typename F>
    void set_returner(const F& function,
                      std::string key = "",
                      bool overwrite = true) {
        if (key.empty()) {
            Processor* p = function();
            key = get_key_and_delete(p);
        }
        ReturnerMap::iterator it = map_.find(key);
        if (it != map_.end()) {
            if (overwrite) {
                it->second = function;
            }
        } else {
            map_[key] = function;
        }
    }

    /** Associate processor type with key.
    \see set_returner()
    */
    template<typename P>
    void set_processor(const std::string& key = "",
                       bool overwrite = true) {
        set_returner(&Meta::new_processor<P>, key, overwrite);
    }

    /** Return keys list */
    Strings keys() const;

    /** Return if no processor returners were set */
    bool empty() const;

    /** Remore all processor returners and options */
    void clear();

    /** Return empty processor which lives till meta object lives.
    This can be used to handle blocksets across different "run"s in script.
    */
    Processor* placeholder_processor() const {
        return placeholder_processor_;
    }

    /** Reset placeholder processor with new instance.
    This can be used to unlink previous actions
    (blocksets) from new ones.
    */
    void reset_placeholder_processor();

    /** Function returning AnyAs */
    typedef boost::function<AnyAs()> AnyReturner;

    /** Get global option */
    AnyAs get_opt(const std::string& key, const AnyAs& dflt = 0) const;

    /** Get option description */
    const std::string& get_description(const std::string& key,
                                       const std::string& dflt = "") const;

    /** Set option description */
    void set_description(const std::string& key,
                         const std::string& description);

    /** Set global option.
    If argument description is empty,
    then description is not changed.
    */
    void set_opt(const std::string& key, const AnyAs& value,
                 const std::string& description = "");

    /** Set global option getter */
    void set_opt_func(const std::string& key, const AnyReturner& f);

    /** Return if global option with this key exists */
    bool has_opt(const std::string& key) const;

    /** List global options */
    Strings opts() const;

    /** Remove global option */
    void remove_opt(const std::string& key);

    lua_State* L() const;

private:
    typedef Processor* ProcessorPtr;
    typedef boost::function<ProcessorPtr()> ProcessorReturner;
    typedef std::map<std::string, ProcessorReturner> ReturnerMap;
    struct GlobalOption {
        AnyReturner f;
        std::string description;
    };
    typedef std::map<std::string, GlobalOption> AnyMap;

    boost::shared_ptr<lua_State> l_;
    // L is initialized before other members
    // L is deleted after other members
    ReturnerMap map_;
    AnyMap opts_;
    Processor* placeholder_processor_;

    template<typename P>
    static P* new_processor() {
        return new P;
    }

    static std::string get_key_and_delete(const Processor* p);
};

}

#endif

