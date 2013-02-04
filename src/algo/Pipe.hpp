/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PIPE_HPP_
#define BR_PIPE_HPP_

#include <vector>

#include "Processor.hpp"

namespace bloomrepeats {

/** Apply several processors */
class Pipe : public Processor {
public:
    /** Constructor */
    Pipe();

    /** Add processor and return *this.
    \note Parameters of processor (workers(), block_set() and other)
        may be changed in this class.
    If block_set is provided, it is set to the processor.
    In run_impl(), if block_set is not set, then block_set of Pipe is set.
    */
    Pipe& add(const ProcessorPtr& processor,
              BlockSetPtr block_set = BlockSetPtr());

    /** Add processor and return *this.
    Overloaded method.
    Ownership is transferred.
    */
    Pipe& add(Processor* processor, BlockSetPtr block_set = BlockSetPtr());

    /** Return max number of applications of all processors.
    All processors are applied in sequence of addition.
    This is repeated untill block set will remain unchanged
    of max_loops() is exceeded.
    -1 means no limit.
    Defaults to 1.
    */
    int max_loops() const {
        return max_loops_;
    }

    /** Set max number of applications of all processors */
    void set_max_loops(int max_loops) {
        max_loops_ = max_loops;
    }

protected:
    /** Add options of all added processors */
    void add_options_impl(po::options_description& desc) const;

    /** Add options to all added processors */
    void apply_options_impl(const po::variables_map& vm);

    /** Apply the action */
    bool run_impl() const;

private:
    std::vector<ProcessorPtr> processors_;
    int max_loops_;
};

}

#endif

