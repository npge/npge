/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_FILE_READER_HPP_
#define BR_FILE_READER_HPP_

#include <iterator>
#include <iosfwd>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace bloomrepeats {

/** Tag for iterator */
class frci_tag;

/** Base class for file readers.
Example:
\code
FileReader reader;
BOOST_FOREACH (std::istream& in, reader) {
    // do smth with in
}

std::istream& in = reader.input(); // get first stream
\endcode
*/
class FileReader {
public:
    /** Iterator class manages file open/close */
    class const_iterator : public std::iterator<frci_tag, std::istream> {
    public:
        /** Destructor */
        ~const_iterator();

        /** Go to next element */
        const_iterator& operator++();

        /** Go to next element */
        const_iterator& operator++(int);

        /** Comparison operator */
        bool operator==(const const_iterator& other);

        /** Comparison operator */
        bool operator!=(const const_iterator& other);

        /** Dereference */
        std::istream& operator*();

    private:
        const FileReader* reader_;
        int index_;
        boost::shared_ptr<std::istream> stream_;

        const_iterator(const FileReader* reader, int index);

        friend class FileReader;
    };

    /** Destructor */
    virtual ~FileReader();

    /** Iterator */
    const_iterator begin() const;

    /** Iterator */
    const_iterator end() const;

    /** Return if no files */
    bool empty() const;

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

    /** First stream.
    If no files, throws Exception.
    */
    std::istream& input() const;

private:
    std::vector<std::string> input_files_;
    mutable boost::shared_ptr<std::istream> stream_;
};

}

#endif

