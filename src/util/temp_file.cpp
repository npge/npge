/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <fstream>
#include <boost/version.hpp>
#if BOOST_VERSION >= 104400
#define BOOST_FILESYSTEM_VERSION 3
#endif
#include <boost/filesystem.hpp>

#include "temp_file.hpp"
#include "Exception.hpp"
#include "reentrant_getenv.hpp"
#include "name_to_stream.hpp"
#include "rand_name.hpp"

namespace npge {

static std::string temp_dir() {
#if BOOST_FILESYSTEM_VERSION == 3
    using namespace boost::filesystem;
    return temp_directory_path().string();
#else
    std::string temp = reentrant_getenv("TEMP");
    if (!temp.empty()) {
        return temp;
    }
    std::string tmpdir = reentrant_getenv("TMPDIR");
    if (!tmpdir.empty()) {
        return tmpdir;
    }
#ifdef _WIN32
    return get_home_dir();
#else
    return "/tmp";
#endif
#endif
}

std::string temp_file() {
    using namespace boost::filesystem;
    using namespace std;
    std::string dir = temp_dir();
#if !defined(_WIN32) && BOOST_FILESYSTEM_VERSION == 3
    const char* const model = "npge-%%%%-%%%%-%%%%-%%%%";
    return unique_path(path(dir) / model).string();
#else
    string result;
    for (int attempt = 0; attempt < 10; attempt++) {
        path p = path(dir) / ("npge_" + rand_name(10));
        string p_s = p.string();
        ofstream file_out(p_s.c_str());
        if (file_out.is_open()) {
            std::string secret = rand_name(10);
            file_out << secret << endl;
            file_out.close();
            if (exists(p_s)) {
                ifstream file_in(p_s.c_str());
                if (file_in.is_open()) {
                    std::string test;
                    file_in >> test;
                    file_in.close();
                    file_out.open(p_s.c_str(),
                                  ios::out | ios::trunc);
                    file_out.close();
                    if (test == secret) {
                        result = p_s;
                        break;
                    }
                }
            }
        }
    }
    return result;
#endif
}

}

