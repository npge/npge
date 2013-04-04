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
#include <boost/assign/list_of.hpp>

#include "name_to_stream.hpp"

namespace bloomrepeats {

// TODO template

using namespace boost::assign; // map_list_of

typedef boost::shared_ptr<std::istream> IstreamPtr;
typedef boost::shared_ptr<std::ostream> OstreamPtr;

static void do_nothing(std::ios_base*)
{ }

static IstreamPtr cin_ptr((&std::cin), do_nothing);
static OstreamPtr cout_ptr((&std::cout), do_nothing);
static OstreamPtr cerr_ptr((&std::cerr), do_nothing);

typedef std::map<std::string, IstreamPtr> Imap;
static Imap custom_istreams_ = map_list_of("", cin_ptr)(":cin", cin_ptr);
static Imap istreams_;
static boost::mutex istreams_mutex_;

typedef std::map<std::string, OstreamPtr> Omap;
static Omap custom_ostreams_ =
    map_list_of("", cout_ptr)(":cout", cout_ptr)(":cerr", cerr_ptr);
static Omap ostreams_;
static boost::mutex ostreams_mutex_;

IstreamPtr name_to_istream(const std::string& name) {
    boost::mutex::scoped_lock lock(istreams_mutex_);
    IstreamPtr stream = istreams_[name];
    if (!stream) {
        Imap::const_iterator it = custom_istreams_.find(name);
        if (it != custom_istreams_.end()) {
            stream = it->second;
        } else if (name.empty() || name[0] == ':') {
            stream.reset(new std::istringstream);
        } else {
            stream.reset(new std::ifstream(name.c_str()));
        }
        istreams_[name] = stream;
    }
    return stream;
}

void set_istream(const std::string& name, IstreamPtr stream) {
    boost::mutex::scoped_lock lock(istreams_mutex_);
    if (stream) {
        custom_istreams_[name] = stream;
    } else {
        custom_istreams_.erase(name);
    }
}

OstreamPtr name_to_ostream(const std::string& name) {
    boost::mutex::scoped_lock lock(ostreams_mutex_);
    OstreamPtr stream = ostreams_[name];
    if (!stream) {
        Omap::const_iterator it = custom_ostreams_.find(name);
        if (it != custom_ostreams_.end()) {
            stream = it->second;
        } else if (name.empty() || name[0] == ':') {
            stream.reset(new std::ostringstream);
        } else {
            stream.reset(new std::ofstream(name.c_str()));
        }
        ostreams_[name] = stream;
    }
    return stream;
}

void set_ostream(const std::string& name, OstreamPtr stream) {
    boost::mutex::scoped_lock lock(ostreams_mutex_);
    if (stream) {
        custom_ostreams_[name] = stream;
    } else {
        custom_ostreams_.erase(name);
    }
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

