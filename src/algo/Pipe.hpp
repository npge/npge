/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PIPE_HPP_
#define BR_PIPE_HPP_

#include "Processor.hpp"

namespace bloomrepeats {

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
    This is repeated untill block set will remain unchanged
    of max_loops() is exceeded.
    -1 means no limit.
    Defaults to 1.
    */
    int max_loops() const;

    /** Set max number of applications of all processors */
    void set_max_loops(int max_loops);

protected:
    /** Add options to all added processors */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

private:
    struct Impl;

    Impl* impl_;
};

}

#endif

