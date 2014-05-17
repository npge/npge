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
#include "boost-xtime.hpp"
#include <boost/thread/mutex.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>

#include "name_to_stream.hpp"

namespace bloomrepeats {

// TODO template

using namespace boost::assign; // map_list_of

typedef boost::shared_ptr<std::istream> IstreamPtr;
typedef boost::shared_ptr<std::ostream> OstreamPtr;

boost::iostreams::stream<boost::iostreams::null_sink> null_ostream(
    (boost::iostreams::null_sink()));

static void do_nothing(std::ios_base*) {
}

static IstreamPtr cin_ptr((&std::cin), do_nothing);
static OstreamPtr cout_ptr((&std::cout), do_nothing);
static OstreamPtr cerr_ptr((&std::cerr), do_nothing);
static OstreamPtr null_ptr((&null_ostream), do_nothing);

typedef std::map<std::string, IstreamPtr> Imap;
static Imap custom_istreams_ = map_list_of("", cin_ptr)(":cin", cin_ptr);
static boost::mutex istreams_mutex_;

typedef std::map<std::string, OstreamPtr> Omap;
static Omap custom_ostreams_ =
    map_list_of("", cout_ptr)(":cout", cout_ptr)
    (":cerr", cerr_ptr)(":null", null_ptr);
static boost::mutex ostreams_mutex_;

IstreamPtr name_to_istream(const std::string& name) {
    boost::mutex::scoped_lock lock(istreams_mutex_);
    Imap::const_iterator it = custom_istreams_.find(name);
    if (it != custom_istreams_.end()) {
        return it->second;
    } else if (name.empty() || name[0] == ':') {
        return IstreamPtr(new std::istringstream);
    } else {
        return IstreamPtr(new std::ifstream(name.c_str()));
    }
}

void set_istream(const std::string& name, IstreamPtr stream) {
    boost::mutex::scoped_lock lock(istreams_mutex_);
    custom_istreams_[name] = stream;
}

void remove_istream(const std::string& name) {
    boost::mutex::scoped_lock lock(istreams_mutex_);
    custom_istreams_.erase(name);
}

OstreamPtr name_to_ostream(const std::string& name) {
    boost::mutex::scoped_lock lock(ostreams_mutex_);
    Omap::const_iterator it = custom_ostreams_.find(name);
    if (it != custom_ostreams_.end()) {
        return it->second;
    } else if (name.empty() || name[0] == ':') {
        return OstreamPtr(new std::ostringstream);
    } else {
        return OstreamPtr(new std::ofstream(name.c_str()));
    }
}

void set_ostream(const std::string& name, OstreamPtr stream) {
    boost::mutex::scoped_lock lock(ostreams_mutex_);
    custom_ostreams_[name] = stream;
}

void remove_ostream(const std::string& name) {
    boost::mutex::scoped_lock lock(ostreams_mutex_);
    custom_ostreams_.erase(name);
}

void set_sstream(const std::string& name, const std::string& c) {
    boost::shared_ptr<std::stringstream> stream(new std::stringstream(c));
    set_istream(name, stream);
    set_ostream(name, stream);
}

void remove_stream(const std::string& name) {
    remove_istream(name);
    remove_ostream(name);
}

void remove_file(const std::string& name) {
    if (!name.empty()) {
        remove(name.c_str());
    }
}

}

