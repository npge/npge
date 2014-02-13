/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FILE_WRITER_HPP_
#define BR_FILE_WRITER_HPP_

#include <string>
#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "global.hpp"

namespace bloomrepeats {

/** Base class for file writers
Empty filename means std::cout.
*/
class FileWriter {
public:
    /** Constructor.
    Call set_remove_after(true)
    */
    FileWriter(Processor* processor, const std::string& opt,
               const std::string& descr, bool required = false);

    /** Destructor.
    Removes the file if get_remove_after().
    */
    virtual ~FileWriter();

    /** Get name of output file */
    std::string output_file() const;

    /** Set name of output file */
    void set_output_file(const std::string& output_file,
                         bool remove_prev = true);

    /** Set random filename inside directory for temporary files */
    void set_rand_name(bool remove_prev = true);

    /** Set if the file will be removed from the destructor */
    void set_remove_after(bool value = true);

    /** Get if the file will be removed from the destructor */
    bool get_remove_after() const;

    /** Return output stream */
    std::ostream& output() const;

private:
    Processor* processor_;
    std::string opt_;
    mutable boost::shared_ptr<std::ostream> output_;
};

}

#endif

