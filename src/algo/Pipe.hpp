/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PIPE_HPP_
#define NPGE_PIPE_HPP_

#include "Processor.hpp"

namespace npge {

/** Apply several processors */
class Pipe : public Processor {
public:
    /** Constructor */
    Pipe(BlockSetPtr other = BlockSetPtr());

    /** Desctructor */
    ~Pipe();

    /** Add processor and return *this.
    Pipe is set as parent of added processor.
    \note Parameters of processor (workers(), block_set() and other)
        may be changed in this class.
    Options is a space separated list of mappings for Processor::set_options().
    */
    Pipe& add(Processor* processor, const std::string& options = "");

    /** Return max number of applications of all processors.
    All processors are applied in sequence of addition.
    This is repeated untill blockset will remain unchanged
    of max_iterations() is exceeded.
    -1 means no limit (until one of previous hashes repeats).
    Defaults to 1.
    */
    int max_iterations() const;

    /** Set max number of applications of all processors */
    void set_max_iterations(int max_iterations);

    /** Return list of processors added */
    std::vector<Processor*> processors() const;

    /** Stop after current iteration */
    void stop();

protected:
    void run_impl() const;

private:
    struct Impl;

    Impl* impl_;
};

}

#endif

