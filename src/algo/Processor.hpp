/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_PROCESSOR_HPP_
#define NPGE_PROCESSOR_HPP_

#include <iosfwd>
#include <string>
#include <boost/utility.hpp>

#include "global.hpp"
#include "po.hpp"
#include "AnyAs.hpp"

namespace npge {

/** Utility class measuring time consumed by processor.
Time-consuming methods of processors must
create instance of this class of stack.
If multiple instances of TimeIncrementer
take same processor instance simultaneouly,
consumed time is calculated as time of last ~TimeIncrementer
minus time of first TimeIncrementer.
*/
class TimeIncrementer {
public:
    /** Constructor */
    TimeIncrementer(const Processor* p);

    /** Destructor */
    ~TimeIncrementer();

private:
    const Processor* p_;
};

/** Wrapper for manipulations with block set */
class Processor : boost::noncopyable {
public:
    /** Function checking options.
    Return value: true if check was passed.
    String: error message or warning.
    */
    typedef boost::function<bool(std::string&)> OptionsChecker;

    /** Function validating option value.
    Gets value of an option and return fixed value.
    */
    typedef boost::function<AnyAs(const AnyAs&)> OptionValidator;

    /** Function returning option value */
    typedef boost::function<AnyAs()> OptionGetter;

    /** Constructor */
    Processor();

    /** Destructor */
    virtual ~Processor();

    /** Declare block set */
    void declare_bs(const std::string& name, const std::string& description);

    /** Remove block set */
    void remove_bs(const std::string& name);

    /** Get description of block set */
    std::string bs_description(const std::string& name) const;

    /** Get named block set.
    If block set with this name is not available,
    create empty block, set it to this name and return.
    */
    BlockSetPtr get_bs(const std::string& name) const;

    /** Set named block set */
    void set_bs(const std::string& name, BlockSetPtr bs);

    /** Return if processor has the block set */
    bool has_bs(const std::string& name) const;

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
        It is usually a Pipe. If processor==0, then parent processor is used
        if defined. If parent processor is not defined and processor==0,
        then blockset rules are not applied.

    Here can be options of several types, separated by space-like chars:
     - mappings, "target=other", see point_bs;
     - default values of properties,
        throug commans line optiond,
        "--workers=10";
        $OPT for global options;
     - ignored option with value, "--workers:=10",
     - "--timing".
     - "no_options".
     - "prefix|prefix-value" see Prefix.

    Note that no spaces are allowed before and after '='.

    Blocks 'target' and 'other' are mapped by default.
    Passing 'target=123' actually means 'target=123 other=other'.
    */
    void set_options(const std::string& options, Processor* processor = 0);

    /** Get "target" block set.
    \deprecated Use get_bs("target").
    */
    BlockSetPtr block_set() const;

    /** Set "target" block set.
    \deprecated Use set_bs("target", block_set).
    */
    void set_block_set(BlockSetPtr block_set);

    /** Get "other" block set.
    \deprecated Use get_bs("other").
    */
    BlockSetPtr other() const;

    /** Set "other" block set.
    \deprecated Use set_bs("other", block_set).
    */
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
    void get_block_sets(Strings& block_sets) const;

    /** Return max number of threads.
    \deprecated Use opt_value("workers").
    */
    int workers() const;

    /** Set max number of threads used to find anchors.
    -1 = number of cores.
    Defaults to 1.
    \deprecated Use set_opt_value("workers", workers).
    */
    void set_workers(int workers);

    /** Add line to log.
    Global variable LOG_TO.
    */
    void write_log(const std::string& message) const;

    /** Close log file */
    void close_log() const;

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
    \deprecated Use opt_value("timing").
    */
    bool timing() const;

    /** Set if this processor prints spent time to stderr from destructor.
    \deprecated Use set_opt_value("timing", timing).
    */
    void set_timing(bool timing);

    /** Copy target and other bs, workers and timing from other processor */
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

    /** Add custom check being run by options_errors() */
    void add_opt_check(const OptionsChecker& checker);

    /** Add custom check in string form.
    \param rule Rule. Syntax: "opt-name operator value-or-other-opt-name".
        Operators: <, >, <=, >=.
        Option types: int, Decimal.
    \param message Error message in case of the check was not passed.
    */
    void add_opt_rule(const std::string& rule, const std::string& message);

    /** Add custom check in string form.
    Shortcut for add_opt_rule(rule, rule).
    */
    void add_opt_rule(const std::string& rule);

    /** Return list of errors with options */
    Strings options_errors() const;

    /** Return list of warnings with options */
    Strings options_warnings() const;

    /** Apply options from strings vector.
    \param options Command line like options.
        Example: ["--workers", "2", "--distance=1"].

    This includes calls to add_options() and apply_options().
    */
    void apply_vector_options(const Strings& options);

    /** Apply options from string.
    \param options Command line like options.
        Example: "--workers 2 --distance=1".

    This includes calls to add_options() and apply_options().
    */
    void apply_string_options(const std::string& options);

    /** Apply the action to the block_set().
    This method calls run_impl() if workers() != 0 && block_set().
    */
    void run() const;

    /** Apply an action to a block.
    Implementation apply_to_block_impl().
    */
    void apply_to_block(Block* block) const;

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
    */
    void apply(const BlockSetPtr& block_set) const;

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

    /** Return clone of this processor.
    New processor is generated by meta with processor key.
    Mappings of blocksets, options, name, parent are copied.
    */
    Processor* clone() const;

    /** Return Meta instance.
    Return meta object set to this processor.
    If it was not set, return meta object of parent.
    If no parent, return 0.

    Meta instance can be set using set_meta().
    */
    Meta* meta() const;

    /** Set Meta object (may be 0).
    Ownership is not transferred.
    */
    void set_meta(Meta* meta);

    /** Get global option value from Meta */
    AnyAs go(const std::string& key, const AnyAs& dflt = 0) const;

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

    /** Add new option to this processor.
    Default value is used to detect type of option.
    Accepted types: int, bool, Decimal, std::string, std::vector<std::string>.
    If type is std::string or std::vector<std::string> and required=true,
    then value of option, provided by user, must not be empty.
    For std::vector<std::string>, default value is ignored.

    If option with such name exists, it is overwritten.
    */
    void add_opt(const std::string& name,
                 const std::string& description,
                 const AnyAs& default_value,
                 bool required = false);

    /** Add new option to this processor.
    This variant sets getter returning value from meta().get_opt(key).
    Option name goes without '$'.
    */
    void add_gopt(const std::string& name,
                  const std::string& description,
                  const std::string& global_opt_name,
                  bool required = false);

    /** Remove option by name.
    \param apply_prefix Whether to apply prefixing
        (add prefix of this processor and its ancestors).
    */
    void remove_opt(const std::string& name, bool apply_prefix = false);

    /** Add option validator.
    Validators are applied to the option value by set_opt_value().
    Requirements:
     - validator must not change default value of the option;
     - result of validator must be of same type as argument;
     - validator must not throw exceptions.
    */
    void add_opt_validator(const std::string& name,
                           const OptionValidator& validator);

    /** Return list of options */
    Strings opts() const;

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
    AnyAs opt_value(const std::string& name) const;

    /** Set value of option.
    \param name Name of option.
    \param value New value of option ($OPT for global option).
    If no option with such name exists, Exception is thrown.
    If type of value differs from type of default value of the option,
    Exception is thrown.
    */
    void set_opt_value(const std::string& name, const AnyAs& value);

    /** Set getter for the option.
    \param name Name of option.
    \param getter New getter of option.
    If getter is set, it is used instead of default value.
    This method does not call the getetr and does not
    check type of value returned by the getter.
    */
    void set_opt_getter(const std::string& name,
                        const OptionGetter& getter);

    /** Set value of option and add it to ignored options */
    void fix_opt_value(const std::string& name, const AnyAs& value);

    /** Set getter for the option and add it to ignored options */
    void fix_opt_getter(const std::string& name,
                        const OptionGetter& getter);

    /** Mark current processor as interrupted.
    When this processor or some of its ruiing children
    notices that processor is_interrupted(), it
    markes the processor and non-interrupted and throws an Exception.
    */
    void interrupt();

    /** Return if this processor or any of its ancestors is interrupted */
    bool is_interrupted() const;

    /** Excape backslashes in input string */
    static std::string escape_backslash(const std::string& str);

    /** Return unique temp file name, removed in destructor.
    \see temp_file()
    */
    std::string tmp_file() const;

protected:
    /** Add options to options description.
    Default implementation does nothing.
    \deprecated Add options in constructor of processor.
    */
    virtual void add_options_impl(po::options_description& desc) const;

    /** Apply options from variables map.
    Default implementation does nothing.
    Implementation may throw Exception.
    \deprecated No need of this. Use opt_value() to get values of options.
    */
    virtual void apply_options_impl(const po::variables_map& vm);

    /** Apply the action to the block_set().
    Default implementation does nothing.
    */
    virtual void run_impl() const;

    /** Apply an action to a block (implementation).
    Creates BlockSet of one block (passed as argument) and call run(block_set).
    */
    virtual void apply_to_block_impl(Block* block) const;

    /** Return human-readable name for the processor.
    Default implementation returns empty line.
    */
    virtual const char* name_impl() const;

    /** Check and process interruption.
    When this processor or some of its ruiing children
    notices that processor is_interrupted(), it
    markes the processor and non-interrupted and throws an Exception.
    */
    void check_interruption() const;

private:
    struct Impl;
    friend class TimeIncrementer;

    Impl* impl_;

    void log_processor(std::ostream& o, int depth);
    void copy_not_ignored(const po::options_description& source,
                          po::options_description& dest) const;
};

/** Return class name by given pointer to processor */
std::string processor_name(const Processor* processor);

/** Streaming operator */
std::ostream& operator<<(std::ostream& o,
                         const Processor& p);

}

#endif

