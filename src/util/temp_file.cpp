/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <boost/version.hpp>
#if BOOST_VERSION >= 104400
#define BOOST_FILESYSTEM_VERSION 3
#endif
#include <boost/filesystem.hpp>

#include "temp_file.hpp"
#include "Exception.hpp"

namespace npge {

static struct Srander {
    Srander() {
        std::srand(time(NULL));
    }
} srander;

std::string temp_file() {
    using namespace boost::filesystem;
    using namespace std;
#if BOOST_FILESYSTEM_VERSION == 3
    const char* const model = "npge-%%%%-%%%%-%%%%-%%%%";
    return unique_path(temp_directory_path() / model).string();
#else
    string result;
    for (int attempt = 0; attempt < 10; attempt++) {
        char file_template[L_tmpnam];
        const char* path_c = tmpnam(file_template);
        if (!path_c) {
            throw Exception("Failed to create temporary file");
        }
        string path(path_c);
        ofstream file_out(path.c_str());
        if (file_out.is_open()) {
            int secret = rand();
            file_out << secret << endl;
            file_out.close();
            if (exists(path)) {
                ifstream file_in(path.c_str());
                if (file_in.is_open()) {
                    int test;
                    file_in >> test;
                    file_in.close();
                    file_out.open(path.c_str(), ios::out | ios::trunc);
                    file_out.close();
                    if (test == secret) {
                        result = path;
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

