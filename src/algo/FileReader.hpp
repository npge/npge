/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_FILE_READER_HPP_
#define NPGE_FILE_READER_HPP_

#include <iterator>
#include <iosfwd>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "global.hpp"

namespace npge {

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
		typedef std::istream value_type;
		typedef int difference_type;
		typedef std::istream* pointer;
		typedef std::istream& reference;
		typedef std::forward_iterator_tag iterator_category;

    public:
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

    /** Constructor */
    FileReader(Processor* processor, const std::string& opt,
               const std::string& descr);

    /** Iterator */
    const_iterator begin() const;

    /** Iterator */
    const_iterator end() const;

    /** Return if no files */
    bool empty() const;

    /** Files list */
    typedef Strings Files;

    /** Get files list */
    Files input_files() const;

    /** Set files list */
    void set_input_files(const Files& input_files);

    /** Set file (list of one file) */
    void set_input_file(const std::string& input_file);

    /** First stream.
    If no files, throws Exception.
    */
    std::istream& input() const;

private:
    Processor* processor_;
    std::string opt_;
    mutable boost::shared_ptr<std::istream> stream_;
};

}

#endif

