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

    /** Destructor */
    virtual ~Processor();

    /** Get target block set */
    const BlockSetPtr& block_set() const {
        return block_set_;
    }

    /** Set target block set */
    void set_block_set(BlockSetPtr block_set) {
        block_set_ = block_set;
    }

    /** Set empty block set */
    void set_empty_block_set();

    /** Return max number of threads */
    int workers() const {
        return workers_;
    }

    /** Set max number of threads used to find anchors.
    -1 = number of cores.
    Defaults to 1.
    */
    void set_workers(int workers);

    /** Get if this processor manages options.
    Defaults to false.
    */
    bool no_options() const {
        return no_options_;
    }

    /** Set if this processor manages options */
    void set_no_options(bool no_options) {
        no_options_ = no_options;
    }

    /** Copy block_set and workers from other processor */
    void assign(const Processor& other);

    /** Add options to options description.
    This method adds --workers option if it was not added yet
    and calls add_options_impl().

    Does nothing if no_options().
    */
    void add_options(po::options_description& desc) const;

    /** Apply options from variables map.
    This method calls apply_options_impl() and
    reads --workers option.
    Implementation may throw Exception.

    Does nothing if no_options().
    */
    void apply_options(const po::variables_map& vm);

    /** Apply the action to the block_set().
    This method calls run_impl() if workers() != 0 && block_set().
    Return if the block set was changed.
    */
    bool run() const;

    /** Return human-readable name for the processor.
    Implementation is name_impl().
    */
    const char* name() const;

    /** Apply the action to other block set.
    This is an equivalent to set_block_set(), run() and set_block_set(previous).
    Return if the block set was changed.
    */
    bool apply(const BlockSetPtr& block_set) const;

protected:
    /** Add options to options description.
    Default implementation does nothing.
    */
    virtual void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map.
    Default implementation does nothing.
    Implementation may throw Exception.
    */
    virtual void apply_options_impl(const po::variables_map& vm);

    /** Apply the action to the block_set().
    Return if the block set was changed.
    Default implementation does nothing.
    */
    virtual bool run_impl() const;

    /** Return human-readable name for the processor.
    Default implementation returns empty line.
    */
    virtual const char* name_impl() const;

private:
    BlockSetPtr block_set_;
    int workers_;
    bool no_options_;

    bool recursive_options() const;
};

}

#endif

