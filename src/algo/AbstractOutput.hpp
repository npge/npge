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

    /** Get output file with all blocks */
    std::string file() const;

    /** Set output file with all blocks */
    void set_file(const std::string& file);

    /** Get mask of output files (${block} is replaced with block name) */
    std::string mask() const;

    /** Set mask of output files */
    void set_mask(const std::string& mask);

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

private:
    std::string file_;
    std::string mask_;
};

}

#endif

