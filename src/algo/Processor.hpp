/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PROCESSOR_HPP_
#define BR_PROCESSOR_HPP_

#include <string>

#include "global.hpp"
#include "po.hpp"

namespace bloomrepeats {

class BlockSetMap;

/** Wrapper for manipulations with block set */
class Processor {
public:
    /** Constructor */
    Processor();

    /** Destructor */
    virtual ~Processor();

    /** Get named block set.
    If block set with this name is not available,
    create empty block, set it to this name and return.
    */
    BlockSetPtr get_bs(const std::string& name) const;

    /** Set named block set */
    void set_bs(const std::string& name, BlockSetPtr bs);

    /** Point named block set to a named block set of other processor.
    Mapping is a string like "target=other" or "target=target",
    where first name is name of block set of this processor,
    and second - of another.
    */
    void point_bs(const std::string& mapping, Processor* processor);

    /** Get "target" block set */
    BlockSetPtr block_set() const;

    /** Set "target" block set */
    void set_block_set(BlockSetPtr block_set);

    /** Get "other" block set */
    BlockSetPtr other() const;

    /** Set "other" block set */
    void set_other(BlockSetPtr other);

    /** Set empty "target" block set */
    void set_empty_block_set();

    /** Set empty "other" block set */
    void set_empty_other();

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

    /** Add option to list of ignored options.
    Ignored options are excluded from options, produced by add_options_impl().

    \warning Ignored options are not passed to apply_options_impl(),
        so it must check their presence.
    */
    void add_ignored_option(const std::string& option);

    /** Get if this processor prints spent time to stderr from destructor.
    Defaults to false.
    */
    bool timing() const {
        return timing_;
    }

    /** Set if this processor prints spent time to stderr from destructor */
    void set_timing(bool timing) {
        timing_ = timing;
    }

    /** Copy block_set and workers from other processor */
    void assign(const Processor& other);

    /** Add options to options description.
    This method adds options --workers, --timing if they were not added yet
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
    If name was set with set_name(), returns that name.
    If name_impl() is not empty, returns name_impl().
    Otherwise returns key().
    */
    std::string name() const;

    /** Set human-readable name for the processor */
    void set_name(const std::string& name) {
        name_ = name;
    }

    /** Apply the action to other block set.
    This is an equivalent to set_block_set(), run() and set_block_set(previous).
    Return if the block set was changed.
    */
    bool apply(const BlockSetPtr& block_set) const;

    /** Return key identifying this class.
    If not key was set, then class name is returned.
    */
    std::string key() const;

    /** Set key identifying this class.
    */
    void set_key(const std::string& key) {
        key_ = key;
    }

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
    BlockSetMap* map_;
    int workers_;
    bool no_options_;
    bool timing_;
    mutable int milliseconds_;
    std::string name_;
    po::options_description ignored_options_;
    std::string key_;
    bool recursive_options() const;
};

/** Return class name by given pointer to processor */
std::string processor_name(const Processor* processor);

/** Return class name by given pointer to processor */
std::string processor_name(const ProcessorPtr& processor);

}

#endif

