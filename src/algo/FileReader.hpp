/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FILE_READER_HPP_
#define BR_FILE_READER_HPP_

#include <string>
#include <vector>

namespace bloomrepeats {

/** Base class for file readers */
class FileReader {
public:
    /** Files list */
    typedef std::vector<std::string> Files;

    /** Get files list */
    const std::vector<std::string>& input_files() const {
        return input_files_;
    }

    /** Set files list */
    void set_input_files(const std::vector<std::string>& input_files) {
        input_files_ = input_files;
    }

    /** Set file (list of one file) */
    void set_input_file(const std::string& input_file);

private:
    std::vector<std::string> input_files_;
};

}

#endif

