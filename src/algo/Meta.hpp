/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_META_HPP_
#define BR_META_HPP_

#include <string>
#include <map>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "global.hpp"

namespace bloomrepeats {

/** Return processor by key */
class Meta {
public:
    /** Constructor.
    Add returners of all processors from algo/.
    */
    Meta();

    /** Return if a processor is associated with the key */
    bool has(const std::string& key) const;

    /** Return processor instance by key.
    If no processor were associated with the key,
    assertation exception is thrown (in debug mode).
    */
    ProcessorPtr get(const std::string& key) const;

    /** Associate function returning processor.
    Key is taken as function()->key().
    */
    template<typename F>
    void set_returner(const F& function) {
        map_[function()->key()] = function;
    }

    /** Associate processor type with key.
    \see set_returner()
    */
    template<typename P>
    void set_processor() {
        set_returner(&Meta::new_processor<P>);
    }

    /** Return keys list */
    std::vector<std::string> keys() const;

    /** Return if no processor returners were set */
    bool empty() const;

    /** Remore all processor returners */
    void clear();

private:
    typedef boost::function<ProcessorPtr()> ProcessorReturner;
    typedef std::map<std::string, ProcessorReturner> ReturnerMap;

    ReturnerMap map_;

    template<typename P>
    static ProcessorPtr new_processor() {
        return ProcessorPtr(new P);
    }
};

}

#endif

