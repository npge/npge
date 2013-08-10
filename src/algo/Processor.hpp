/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef BR_PROCESSOR_HPP_
#define BR_PROCESSOR_HPP_

#include <iosfwd>
#include <string>
#include <boost/utility.hpp>

#include "global.hpp"
#include "po.hpp"
#include "AnyAs.hpp"

namespace bloomrepeats {

/** Wrapper for manipulations with block set */
class Processor : boost::noncopyable {
public:
    /** Function checking options.
    Return value: true of check was passed.
    String: error message or warning.
    */
    typedef boost::function<bool(std::string&)> OptionsChecker;

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
    \param mapping String like "target=other" or "target=target",
        where first name is name of block set of this processor,
        and second - of another.
    \param processor Processor from which blocks are taken by name.
        It is usually a Pipe.
    */
    void point_bs(const std::string& mapping, Processor* processor);

    /** Set text options.
    \param options String with options.
        where first name is name of block set of this processor,
        and second - of another.
    \param processor Processor from which blocks are taken by name.
        It is usually a Pipe. If processor==0, then blockset mapping
        is not allowed.

    Here can be options of several types, separated by space-like chars:
     - mappings, "target=other", see point_bs;
     - default values of properties, throug commans line optiond,
        "--workers=10";
     - ignored option with value, "--workers:=10",
     - "--timing".
     - "no_options".
     - "no_remove_after" see FileWriter.
     - "prefix|prefix-value" see Prefix.

    Note that no spaces are allowed before and after '='.

    Blocks 'target' and 'other' are mapped by default.
    Passing 'target=123' actually means 'target=123 other=other'.
    */
    void set_options(const std::string& options, Processor* processor = 0);

    /** Get "target" block set */
    BlockSetPtr block_set() const;

    /** Set "target" block set */
    void set_block_set(BlockSetPtr block_set);

    /** Get "other" block set */
    BlockSetPtr other() const;

    /** Set "other" block set */
    void set_other(BlockSetPtr other);

    /** Set empty "target" block set.
    \deprecated Useless method. Empty bloc set is created automagicaly.
    */
    void set_empty_block_set();

    /** Set empty "other" block set.
    \deprecated Useless method. Empty bloc set is created automagicaly.
    */
    void set_empty_other();

    /** Appends names of block sets to vector */
    void get_block_sets(std::vector<std::string>& block_sets) const;

    /** Return max number of threads */
    int workers() const;

    /** Set max number of threads used to find anchors.
    -1 = number of cores.
    Defaults to 1.
    */
    void set_workers(int workers);

    /** Get if this processor manages options.
    Defaults to false.
    */
    bool no_options() const;

    /** Set if this processor manages options */
    void set_no_options(bool no_options);

    /** Add option to list of ignored options.
    Ignored options are excluded from options, produced by add_options_impl().

    \warning Ignored options are not passed to apply_options_impl(),
        so it must check their presence.
    */
    void add_ignored_option(const std::string& option);

    /** Return if the option is ignored by this Processor or its ancestors */
    bool is_ignored(const std::string& option) const;

    /** Get if this processor prints spent time to stderr from destructor.
    Defaults to false.
    */
    bool timing() const;

    /** Set if this processor prints spent time to stderr from destructor */
    void set_timing(bool timing);

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

    Does nothing except --timing and --workers if no_options().
    */
    void apply_options(const po::variables_map& vm);

    /** Return list of errors with options */
    std::vector<std::string> options_errors() const;

    /** Return list of warnings with options */
    std::vector<std::string> options_warnings() const;

    /** Apply options from string.
    \param options Command line like options.
        Example: ["--workers", "2", "--distance=1"].

    This includes calls to add_options() and apply_options().
    */
    void apply_vector_options(const std::vector<std::string>& options);

    /** Apply options from string.
    \param options Command line like options.
        Example: "--workers 2 --distance=1".

    This includes calls to add_options() and apply_options().
    */
    void apply_string_options(const std::string& options);

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
    void set_name(const std::string& name);

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
    void set_key(const std::string& key);

    /** Return parent processor (may be 0) */
    Processor* parent() const;

    /** Set parent processor (may be 0).
    Parent processor owns children.
    */
    void set_parent(Processor* parent);

    /** Return list of children */
    std::vector<Processor*> children() const;

    /** Return Meta instance.
    Return meta object set to this processor.
    If it was not set, return meta object of parent.
    If no parent, return default Meta (static thread-specific).

    Meta instance can be set using set_meta().
    */
    Meta* meta() const;

    /** Set Meta object (may be 0).
    Ownership is not transferred.
    */
    void set_meta(Meta* meta);

    /** Get options prefix of this processor.
    This string is added before names of options of this processor.
    If processor has parent, its prefix added before prefix of this processor,
    parent's parent before all of them, etc.
    Example: "hits-"
    Default prefix is empty string.
    */
    const std::string& opt_prefix() const;

    /** Set options prefix of this processor */
    void set_opt_prefix(const std::string& prefix);

    /** Return options-prefixed string.
    Add prefixes of this processor and all its ancestors
    before the string and return result.
    \note Setting a parent or setting prefix (or prefix of
    some parent) may affect result of this function.
    */
    std::string opt_prefixed(const std::string& name) const;

    /** Return list of options */
    std::vector<std::string> opts() const;

    /** Return if the processor has option with this name.
    \param name Name of option.
    */
    bool has_opt(const std::string& name) const;

    /** Return description of the option.
    \param name Name of option.
    If no option with such name exists, Exception is thrown.
    */
    const std::string& opt_description(const std::string& name) const;

    /** Return type of the option.
    \param name Name of option.
    If no option with such name exists, Exception is thrown.
    */
    const std::type_info& opt_type(const std::string& name) const;

    /** Return default value of option.
    \param name Name of option.
    If no option with such name exists, Exception is thrown.
    */
    const AnyAs& default_opt_value(const std::string& name) const;

    /** Return value of option.
    \param name Name of option.
    If no option with such name exists, Exception is thrown.
    */
    const AnyAs& opt_value(const std::string& name) const;

    /** Set value of option.
    \param name Name of option.
    \param value New value of option.
    If no option with such name exists, Exception is thrown.
    If type of value differs from type of default value of the option,
    Exception is thrown.
    */
    void set_opt_value(const std::string& name, const AnyAs& value);

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

    /** Add new option to this processor.
    Default value is used to detect type of option.
    Accepted types: int, bool, double, std::string, std::vector<std::string>.
    If type is std::string or std::vector<std::string> and required=true,
    then value of option, provided by user, must not be empty.
    For std::vector<std::string>, default value is ignored.

    If option with such name exists, it is overwritten.
    */
    void add_opt(const std::string& name,
                 const std::string& description,
                 const AnyAs& default_value,
                 bool required = false);

    /** Remove option by name.
    \param apply_prefix Whether to apply prefixing
        (add prefix of this processor and its ancestors).
    */
    void remove_opt(const std::string& name, bool apply_prefix = false);

    /** Add custom check being run by options_errors() */
    void add_options_check(const OptionsChecker& checker);

    /** Add custom check in string form.
    \param rule Rule. Syntax: "opt-name operator value-or-other-opt-name".
        Operators: <, >, <=, >=.
        Option types: int, double.
    \param message Error message in case of the check was not passed.
    */
    void add_options_check(const std::string& rule, const std::string& message);

    /** Add custom check in string form.
    Shortcut for add_options_check(rule, rule).
    */
    void add_options_check(const std::string& rule);

private:
    struct Impl;

    Impl* impl_;

    void log_processor(std::ostream& o, int depth);
    void copy_not_ignored(const po::options_description& source,
                          po::options_description& dest) const;
};

/** Return class name by given pointer to processor */
std::string processor_name(const Processor* processor);

}

#endif

