/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#define NPG_UNIX
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#endif

#include <cstdlib>
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

namespace npge {

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

static std::string resolve_home_dir(const std::string& d) {
    if (!d.empty() && d[0] == '~') {
        return get_home_dir() + d.substr(1);
    } else {
        return d;
    }
}

IstreamPtr name_to_istream(const std::string& name) {
    boost::mutex::scoped_lock lock(istreams_mutex_);
    Imap::const_iterator it = custom_istreams_.find(name);
    if (it != custom_istreams_.end()) {
        return it->second;
    } else if (name.empty() || name[0] == ':') {
        return IstreamPtr(new std::istringstream);
    } else {
        std::string path = resolve_home_dir(name);
        return IstreamPtr(new std::ifstream(path.c_str()));
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

// http://stackoverflow.com/a/13062069
template<typename T>
struct array_deleter {
    void operator()(T const* p) {
        delete[] p;
    }
};

std::string get_home_dir(const std::string& dftl) {
    // http://stackoverflow.com/a/3733955
#ifdef NPG_UNIX
    // Unix
    const char* home = getenv("HOME");
    if (home) {
        return home;
    } else {
#if _POSIX_C_SOURCE >= 1 || _XOPEN_SOURCE || _BSD_SOURCE || _SVID_SOURCE || _POSIX_SOURCE
        // getpwuid_r
        int size = sysconf(_SC_GETPW_R_SIZE_MAX);
        if (size == -1) {
            size = 1024;
        }
        boost::shared_ptr<char> buf(new char[size],
                                    array_deleter<char>());
        passwd* p_ptr;
        passwd p;
        getpwuid_r(getuid(), &p, buf.get(), size, &p_ptr);
        if (p_ptr) {
            return p_ptr->pw_dir;
        }
#else // getpwuid
        struct passwd* pw = getpwuid(getuid());
        return pw->pw_dir;
#endif
    }
#else
    // Windows
    const char* userprofile = getenv("USERPROFILE");
    if (userprofile) {
        return userprofile;
    } else {
        const char* homedrive = getenv("HOMEDRIVE");
        const char* homepath = getenv("HOMEPATH");
        if (homedrive && homepath) {
            return std::string(homedrive) + homepath;
        }
    }
#endif
    return dftl;
}

}

