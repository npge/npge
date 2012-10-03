/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PROCESSOR_HPP_
#define BR_PROCESSOR_HPP_

#include "global.hpp"
#include "po.hpp"

namespace bloomrepeats {

/** Wrapper for manipulations with block set */
class Processor {
public:
    /** Constructor */
    Processor();

    /** Get target block set */
    const BlockSetPtr& block_set() const {
        return block_set_;
    }

    /** Set target block set */
    void set_block_set(BlockSetPtr block_set) {
        block_set_ = block_set;
    }

    /** Return max number of threads */
    int workers() const {
        return workers_;
    }

    /** Set max number of threads used to find anchors.
    -1 = number of cores.
    Defaults to 1.
    */
    void set_workers(int workers);

    /** Copy block_set and workers from other processor */
    void assign(const Processor& other);

    /** Add options to options description.
    This method calls add_options_impl() and
    adds --workers option if it was not added yet.
    */
    void add_options(po::options_description& desc) const;

    /** Apply options from variables map.
    This method calls apply_options_impl() and
    reads --workers option.
    */
    void apply_options(const po::variables_map& vm);

    /** Apply the action to the block_set().
    This method calls run_impl() if workers() != 0 && block_set().
    */
    void run();

protected:
    /** Add options to options description.
    Default implementation does nothing.
    */
    virtual void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map.
    Default implementation does nothing.
    */
    virtual void apply_options_impl(const po::variables_map& vm);

    /** Apply the action to the block_set().
    Default implementation does nothing.
    */
    virtual void run_impl();

private:
    BlockSetPtr block_set_;
    int workers_;
};

}

#endif

