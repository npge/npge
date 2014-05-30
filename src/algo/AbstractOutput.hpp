/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_ABSTRACT_OUTPUT_HPP_
#define NPGE_ABSTRACT_OUTPUT_HPP_

#include <iosfwd>

#include "BlocksJobs.hpp"

namespace npge {

/** Print some function from blocks to file or to stdout */
class AbstractOutput : public BlocksJobs {
public:
    /** Constructor */
    AbstractOutput();

    /** Destructor */
    ~AbstractOutput();

    /** Return if output for all blocks is written to one file */
    bool one_file() const;

protected:
    /** Sort blocks by size desc, then by name */
    void change_blocks_impl(std::vector<Block*>& blocks) const;

    void initialize_work_impl() const;

    ThreadData* before_thread_impl() const;

    void process_block_impl(Block* block, ThreadData* data) const;

    void finish_work_impl() const;

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
    struct Impl;
    Impl* impl_;
};

}

#endif

