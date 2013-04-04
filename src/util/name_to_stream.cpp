/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdio>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/thread/mutex.hpp>

#include "name_to_stream.hpp"

namespace bloomrepeats {

// TODO template

typedef boost::shared_ptr<std::istream> IstreamPtr;
typedef boost::shared_ptr<std::ostream> OstreamPtr;

static std::map<std::string, IstreamPtr> istreams_;
static boost::mutex istreams_mutex_;

static std::map<std::string, OstreamPtr> ostreams_;
static boost::mutex ostreams_mutex_;

static void do_nothing(std::ios_base*)
{ }

static IstreamPtr cin_ptr((&std::cin), do_nothing);
static OstreamPtr cout_ptr((&std::cout), do_nothing);

IstreamPtr name_to_istream(const std::string& name) {
    if (name.empty()) {
        return cin_ptr;
    }
    boost::mutex::scoped_lock lock(istreams_mutex_);
    IstreamPtr stream = istreams_[name];
    if (!stream) {
        if (name[0] == ':') {
            stream.reset(new std::istringstream);
        } else {
            stream.reset(new std::ifstream(name.c_str()));
        }
        istreams_[name] = stream;
    }
    return stream;
}

OstreamPtr name_to_ostream(const std::string& name) {
    if (name.empty()) {
        return cout_ptr;
    }
    boost::mutex::scoped_lock lock(ostreams_mutex_);
    OstreamPtr stream = ostreams_[name];
    if (!stream) {
        if (name[0] == ':') {
            stream.reset(new std::ostringstream);
        } else {
            stream.reset(new std::ofstream(name.c_str()));
        }
        ostreams_[name] = stream;
    }
    return stream;
}

void remove_istream(const std::string& name) {
    if (!name.empty()) {
        boost::mutex::scoped_lock lock(istreams_mutex_);
        istreams_.erase(name);
    }
}

void remove_ostream(const std::string& name, bool remove_file) {
    if (!name.empty()) {
        boost::mutex::scoped_lock lock(ostreams_mutex_);
        ostreams_.erase(name);
        if (name[0] != ':' && remove_file) {
            remove(name.c_str());
        }
    }
}

}

