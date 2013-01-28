/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_OPTIONS_PREFIX_HPP
#define BR_OPTIONS_PREFIX_HPP

#include <string>

#include "po.hpp"

namespace bloomrepeats {

/** Adds a prefix to each option of descedant class */
class OptionsPrefix {
public:
    /** Constructor */
    OptionsPrefix(const std::string& prefix = "");

    /** Get prefix */
    const std::string& prefix() const {
        return prefix_;
    }

    /** Get prefixed string */
    std::string prefixed(const std::string& option) const;

    /** Set prefix */
    void set_prefix(const std::string& prefix) {
        prefix_ = prefix;
    }

protected:
    class PrefixOptions : public AddUniqueOptions {
    public:
        PrefixOptions(po::options_description& desc,
                      const std::string* prefix);

        PrefixOptions& operator()(const char* name,
                                  const char* description);

        PrefixOptions& operator()(const char* name,
                                  const po::value_semantic* s,
                                  const char* description);

        PrefixOptions& operator()(const char* name,
                                  const po::value_semantic* s);

    private:
        const std::string* prefix_;
    };

    PrefixOptions add_unique_options(po::options_description& desc) const;

private:
    std::string prefix_;
};

}

#endif

