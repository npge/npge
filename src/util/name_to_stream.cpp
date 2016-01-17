/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#define NPG_UNIX
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#elif defined(_WIN32) || defined(__WIN32__)
#include <winsock2.h>
#include <windows.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include "boost-xtime.hpp"
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>
#if BOOST_VERSION >= 104400
#define BOOST_FILESYSTEM_VERSION 3
#else
#define BOOST_FILESYSTEM2_NARROW_ONLY
#endif
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "name_to_stream.hpp"
#include "reentrant_getenv.hpp"
#include "Exception.hpp"

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
static Imap custom_istreams_ = map_list_of("", cin_ptr)(":stdin", cin_ptr);
static boost::mutex istreams_mutex_;

typedef std::map<std::string, OstreamPtr> Omap;
static Omap custom_ostreams_ =
    map_list_of("", cout_ptr)(":stdout", cout_ptr)
    (":stderr", cerr_ptr)(":null", null_ptr);
static boost::mutex ostreams_mutex_;

std::string resolve_home_dir(const std::string& d0) {
    using namespace boost::algorithm;
    std::string d = replace_first_copy(d0, "{app_dir}",
                                       get_app_dir());
    if (!d.empty() && d[0] == '~') {
        return get_home_dir() + d.substr(1);
    } else {
        return d;
    }
}

// http://stackoverflow.com/a/6821180

const int BUFFER_SIZE = 4096;

struct BufferedIfstream : public std::ifstream {
    char buffer_[BUFFER_SIZE];
};

IstreamPtr name_to_istream(const std::string& name) {
    boost::mutex::scoped_lock lock(istreams_mutex_);
    Imap::const_iterator it = custom_istreams_.find(name);
    if (it != custom_istreams_.end()) {
        return it->second;
    } else if (name.empty() || name[0] == ':') {
        return boost::make_shared<std::istringstream>();
    } else {
        std::string path = resolve_home_dir(name);
        boost::shared_ptr<BufferedIfstream> result =
            boost::make_shared<BufferedIfstream>();
        result->rdbuf()->pubsetbuf(result->buffer_, BUFFER_SIZE);
        result->open(path.c_str(),
            std::ios_base::in | std::ios_base::binary);
        if (!result->is_open()) {
            throw Exception("Error opening file " + name);
        }
        return result;
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
        return boost::make_shared<std::ostringstream>();
    } else {
        boost::shared_ptr<std::ofstream> result =
            boost::make_shared<std::ofstream>(name.c_str());
        if (!result->is_open()) {
            throw Exception("Error opening file " + name);
        }
        return result;
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
    typedef boost::shared_ptr<std::stringstream> SPtr;
    SPtr stream = boost::make_shared<std::stringstream>(c);
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
    std::string home = reentrant_getenv("HOME");
    if (!home.empty()) {
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
    std::string userprofile = reentrant_getenv("USERPROFILE");
    if (!userprofile.empty()) {
        return userprofile;
    } else {
        std::string homedrive = reentrant_getenv("HOMEDRIVE");
        std::string homepath = reentrant_getenv("HOMEPATH");
        if (!homedrive.empty() && !homepath.empty()) {
            return homedrive + homepath;
        }
    }
#endif
    return dftl;
}

namespace fs = boost::filesystem;

std::string get_current_dir(const std::string&) {
    return fs::current_path().string();
}

char* app_path_ = 0;
boost::mutex app_path_mutex_;

struct AppPathDeleter {
    ~AppPathDeleter() {
        delete[] app_path_;
    }
} app_path_deleter_;

std::string get_app_path() {
    boost::mutex::scoped_lock lock(app_path_mutex_);
    if (!app_path_) {
        const int MAX_APP_PATH = 1000;
        char path[MAX_APP_PATH];
        int s = 0;
#ifdef NPG_UNIX
        s = readlink("/proc/self/exe", path, MAX_APP_PATH);
#elif defined(_WIN32) || defined(__WIN32__)
        HMODULE hModule = GetModuleHandle(NULL);
        if (hModule != NULL) {
            s = GetModuleFileName(hModule, path, MAX_APP_PATH);
        }
#endif
        if (s > 0) {
            app_path_ = new char[s + 1];
            memcpy(app_path_, path, s);
            app_path_[s] = '\0';
        } else {
            // UPX preserves content of readlink("/proc/self/exe")
            // in environment variable named " " [three spaces]
            std::string upx_self_exe = reentrant_getenv("   ");
            app_path_ = new char[upx_self_exe.size() + 1];
            strcpy(app_path_, upx_self_exe.c_str());
        }
    }
    return app_path_;
}

std::string get_app_dir() {
    return fs::path(get_app_path()).parent_path().string();
}

std::string cat_paths(const std::string& path1,
                      const std::string& path2) {
    return (fs::path(path1) / fs::path(path2)).string();
}

std::string system_complete(const std::string& p) {
    return fs::system_complete(fs::path(p)).string();
}

bool file_exists(const std::string& p) {
    return fs::exists(fs::path(p));
}

bool is_dir(const std::string& p) {
    return fs::is_directory(fs::path(p));
}

Strings dir_children(const std::string& d) {
    Strings result;
    fs::directory_iterator dir(d), end;
    BOOST_FOREACH (const fs::path& child,
                  std::make_pair(dir, end)) {
        result.push_back(child.string());
    }
    return result;
}

void copy_file(const std::string& src,
               const std::string& dst) {
#if BOOST_FILESYSTEM_VERSION == 3
    fs::copy_file(src, dst,
                  fs::copy_option::overwrite_if_exists);
#else
    fs::copy_file(fs::path(src), fs::path(dst),
                  fs::copy_option::overwrite_if_exists);
#endif
}

void make_dir(const std::string& dir) {
    fs::create_directory(dir);
}

std::string to_filename(const std::string& p) {
#if BOOST_FILESYSTEM_VERSION == 3
    return fs::path(p).filename().string();
#else
    return fs::path(p).filename();
#endif
}

std::string escape_path(const std::string& str) {
    using namespace boost::algorithm;
    return replace_all_copy(str, "\\", "\\\\");
}

}

