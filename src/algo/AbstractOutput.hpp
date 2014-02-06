/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_ABSTRACT_OUTPUT_HPP_
#define BR_ABSTRACT_OUTPUT_HPP_

#include <iosfwd>

#include "Processor.hpp"
#include "OptionsPrefix.hpp"

namespace bloomrepeats {

/** Print some function from blocks to file or to stdout */
class AbstractOutput : public Processor {
public:
    /** Constructor */
    AbstractOutput();

    /** Return if output for all blocks is written to one file */
    bool one_file() const;

protected:
    bool run_impl() const;

    /** Do something in the beginning.
    By default, does nothing.
    */
    virtual void prepare() const;

    /** Print header of output file (once per file).
    By default, does nothing.
    */
    virtual void print_header(std::ostream& o) const;

    /** Print block to output stream */
    virtual void print_block(std::ostream& o, Block* block) const = 0;

    /** Print footer of output file (once per file).
    By default, does nothing.
    */
    virtual void print_footer(std::ostream& o) const;
};

}

#endif

