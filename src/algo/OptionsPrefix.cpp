/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "OptionsPrefix.hpp"

namespace bloomrepeats {

OptionsPrefix::OptionsPrefix(const std::string& prefix):
    prefix_(prefix)
{ }

std::string OptionsPrefix::prefixed(const std::string& option) const {
    return prefix() + option;
}

OptionsPrefix::PrefixOptions::PrefixOptions(po::options_description& desc,
        const std::string* prefix):
    AddUniqueOptions(desc),
    prefix_(prefix)
{ }

typedef po::options_description_easy_init easy_init;
#define PO_ OptionsPrefix::PrefixOptions

PO_& PO_::operator()(const char* name,
                     const char* description) {
    easy_init::operator()((*prefix_ + name).c_str(), description);
    return *this;
}

PO_& PO_::operator()(const char* name,
                     const po::value_semantic* s,
                     const char* description) {
    easy_init::operator()((*prefix_ + name).c_str(), s, description);
    return *this;
}

PO_& PO_::operator()(const char* name,
                     const po::value_semantic* s) {
    easy_init::operator()((*prefix_ + name).c_str(), s);
    return *this;
}

PO_ OptionsPrefix::add_unique_options(po::options_description& desc) const {
    return PO_(desc, &prefix_);
}

}

