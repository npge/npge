/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FILE_WRITER_HPP_
#define BR_FILE_WRITER_HPP_

#include <string>

namespace bloomrepeats {

/** Base class for file writers */
class FileWriter {
public:
    /** Constructor.
    Call set_rand_name() and set_remove_after(true)
    */
    FileWriter();

    /** Destructor.
    Removes the file if get_remove_after().
    */
    virtual ~FileWriter();

    /** Get name of output file */
    const std::string& file() const {
        return file_;
    }

    /** Set name of output file */
    void set_file(const std::string& file, bool remove_prev = true);

    /** Set random filename inside directory for temporary files */
    void set_rand_name(bool remove_prev = true);

    /** Remove file */
    void remove_file();

    /** Set if the file will be removed from the destructor */
    void set_remove_after(bool value = true);

    /** Get if the file will be removed from the destructor */
    bool get_remove_after() const {
        return remove_after_;
    }

private:
    std::string file_;
    bool remove_after_;
};

}

#endif

