/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <typeinfo>
#include "boost-xtime.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "Processor.hpp"
#include "BlockSet.hpp"
#include "FileWriter.hpp"
#include "class_name.hpp"
#include "string_arguments.hpp"
#include "throw_assert.hpp"
#include "tss_meta.hpp"
#include "Exception.hpp"
#include "name_to_stream.hpp"

namespace bloomrepeats {

struct BlockSetHolder {
public:
    BlockSetHolder():
        processor_(0)
    { }

    BlockSetPtr block_set() const {
        if (block_set_) {
            return block_set_;
        } else if (processor_) {
            return processor_->get_bs(name_);
        } else {
            block_set_ = new_bs();
            return block_set_;
        }
    }

    void set_block_set(BlockSetPtr block_set) {
        block_set_ = block_set;
        processor_ = 0;
        name_.clear();
    }

    void set_processor(Processor* processor, const std::string& name) {
        block_set_.reset();
        processor_ = processor;
        name_ = name;
    }

    std::string description_;

private:
    mutable BlockSetPtr block_set_;
    // or
    Processor* processor_;
    std::string name_;
};

typedef std::map<std::string, BlockSetHolder> BlockSetMap;

struct Option {
    Option():
        required_(false)
    { }

    Option(const std::string& name,
           const std::string& description,
           const AnyAs& default_value = AnyAs(),
           bool required = false):
        name_(name),
        description_(description),
        default_value_(default_value),
        required_(required)
    { }

    std::string name_;
    std::string description_;
    AnyAs default_value_;
    AnyAs value_;
    bool required_;
    std::vector<Processor::OptionValidator> validators_;

    const std::type_info& type() const {
        return default_value_.type();
    }

    const AnyAs& final_value() const {
        if (!value_.empty()) {
            return value_;
        } else {
            return default_value_;
        }
    }
};

struct Processor::Impl {
    Impl():
        no_options_(false), milliseconds_(0),
        logged_(false), parent_(0), meta_(0)
    { }

    BlockSetMap map_;
    bool no_options_;
    mutable int milliseconds_;
    std::string name_;
    po::options_description ignored_options_;
    std::string key_;
    bool logged_;
    Processor* parent_;
    std::vector<Processor*> children_;
    Meta* meta_;
    std::string opt_prefix_;
    typedef std::map<std::string, Option> Name2Option; // option name to Option
    Name2Option opts_;
    std::vector<OptionsChecker> checkers_;
};

static AnyAs workers_1(AnyAs workers) {
    int value = workers.as<int>();
    if (value == -1) {
        value = boost::thread::hardware_concurrency();
    }
    return value;
}

Processor::Processor() {
    impl_ = new Impl;
    add_opt("workers", "number of threads", 1);
    add_opt_validator("workers", workers_1);
    add_opt("timing", "measure time for each processor", false);
}

Processor::~Processor() {
    if (!impl_->logged_ && impl_->milliseconds_) {
        log_processor(*name_to_ostream(":cerr"), 0);
    }
    BOOST_FOREACH (Processor* child, impl_->children_) {
        child->impl_->parent_ = 0;
        delete child;
    }
    impl_->children_.clear();
    if (parent()) {
        set_parent(0);
    }
    delete impl_;
}

void Processor::declare_bs(const std::string& name,
                           const std::string& description) {
    impl_->map_[name].description_ = description;
}

void Processor::remove_bs(const std::string& name) {
    impl_->map_.erase(name);
}

std::string Processor::bs_description(const std::string& name) const {
    BlockSetMap::const_iterator it = impl_->map_.find(name);
    if (it == impl_->map_.end()) {
        return "";
    } else {
        return it->second.description_;
    }
}

BlockSetPtr Processor::get_bs(const std::string& name) const {
    BlockSetMap::const_iterator it = impl_->map_.find(name);
    if (it == impl_->map_.end()) {
        BlockSetPtr bs = new_bs();
        impl_->map_[name].set_block_set(bs);
        return bs;
    } else {
        return it->second.block_set();
    }
}

void Processor::set_bs(const std::string& name, BlockSetPtr bs) {
    impl_->map_[name].set_block_set(bs);
}

void Processor::point_bs(const std::string& mapping, Processor* processor) {
    size_t eq_pos = mapping.find("=");
    BOOST_ASSERT_MSG(eq_pos != std::string::npos,
                     ("Bad mapping: " + mapping).c_str());
    std::string name_in_this = mapping.substr(0, eq_pos);
    std::string name_in_processor = mapping.substr(eq_pos + 1);
    BOOST_ASSERT_MSG(processor != this || name_in_this != name_in_processor,
                     ("Trying to set self-pointed blockset: " + mapping +
                      " in processor " + key()).c_str());
    impl_->map_[name_in_this].set_processor(processor, name_in_processor);
}

void Processor::set_options(const std::string& options, Processor* processor) {
    if (!processor) {
        processor = parent();
    }
    if (processor) {
        point_bs("target=target", processor);
        point_bs("other=other", processor);
    }
    std::vector<std::string> ignored;
    std::vector<std::string> default_opts;
    using namespace boost::algorithm;
    using boost::tokenizer;
    using boost::escaped_list_separator;
    typedef tokenizer<escaped_list_separator<char> > tok_t;
    tok_t tok(options, escaped_list_separator<char>('\\', ' ', '\"'));
    BOOST_FOREACH (std::string opt, tok) {
        trim_right(opt);
        size_t eq_pos = opt.find('=');
        if (eq_pos != std::string::npos) {
            if (opt[0] == '-') {
                // command line option
                std::string opt_name = opt.substr(0, eq_pos);
                std::string opt_value = opt.substr(eq_pos + 1);
                bool ignore = opt_name[opt_name.size() - 1] == ':';
                if (ignore) {
                    opt_name.resize(opt_name.size() - 1);
                }
                default_opts.push_back(opt_name);
                default_opts.push_back(opt_value);
                if (ignore) {
                    std::string short_name;
                    if (opt_name.size() > 1 && opt_name[1] == '-') {
                        // option like --workers
                        short_name = opt_name.substr(2, opt_name.size() - 2);
                    } else {
                        // option like -i
                        short_name = opt_name.substr(1, opt_name.size() - 1);
                    }
                    ignored.push_back(short_name);
                }
            } else if (processor) {
                point_bs(opt, processor);
            } else {
                // TODO bad option
            }
        } else if (opt == "no_options") {
            set_no_options(true);
        } else if (starts_with(opt, "prefix|")) {
            int sep = opt.find('|');
            std::string prefix_value = opt.substr(sep + 1);
            set_opt_prefix(prefix_value);
        } else {
            // TODO bad option
        }
    }
    if (!default_opts.empty()) {
        apply_vector_options(default_opts);
    }
    BOOST_FOREACH (const std::string& opt, ignored) {
        add_ignored_option(opt);
    }
}

BlockSetPtr Processor::block_set() const {
    return get_bs("target");
}

void Processor::set_block_set(BlockSetPtr block_set) {
    set_bs("target", block_set);
}

BlockSetPtr Processor::other() const {
    return get_bs("other");
}

void Processor::set_other(BlockSetPtr other) {
    set_bs("other", other);
}

void Processor::set_empty_block_set() {
    set_block_set(new_bs());
}

void Processor::set_empty_other() {
    set_other(new_bs());
}

void Processor::get_block_sets(std::vector<std::string>& block_sets) const {
    BOOST_FOREACH (const BlockSetMap::value_type& name_and_bs, impl_->map_) {
        const std::string& name = name_and_bs.first;
        block_sets.push_back(name);
    }
}

int Processor::workers() const {
    return opt_value("workers").as<int>();
}

void Processor::set_workers(int workers) {
    set_opt_value("workers", workers);
}

bool Processor::no_options() const {
    return impl_->no_options_;
}

void Processor::set_no_options(bool no_options) {
    impl_->no_options_ = no_options;
}

void Processor::add_ignored_option(const std::string& option) {
    add_unique_options(impl_->ignored_options_)(option.c_str(), "");
}

bool Processor::is_ignored(const std::string& option) const {
    const Processor* p = this;
    std::string name = option;
    while (p) {
        name = p->opt_prefix() + name;
        if (p->impl_->ignored_options_.find_nothrow(name, /* approx */ false)) {
            return true;
        }
        p = p->parent();
    }
    return false;
}

bool Processor::timing() const {
    return opt_value("timing").as<bool>();
}

void Processor::set_timing(bool timing) {
    set_opt_value("timing", timing);
}

void Processor::assign(const Processor& other) {
    set_block_set(other.block_set());
    set_other(other.other());
    set_workers(other.workers());
    set_timing(other.timing());
}

static bool good_opt_type(const std::type_info& ti) {
    return ti == typeid(int) || ti == typeid(bool) || ti == typeid(double) ||
           ti == typeid(std::string) || ti == typeid(std::vector<std::string>);
}

static void add_option(po::options_description& desc, const std::string name,
                       const Option& opt) {
    BOOST_ASSERT(good_opt_type(opt.type()));
    typedef boost::shared_ptr<po::option_description> OptPtr;
    po::value_semantic* vs = 0;
    if (opt.type() == typeid(int)) {
        vs = po::value<int>()->default_value(opt.final_value().as<int>());
    } else if (opt.type() == typeid(bool)) {
        vs = po::value<bool>()->default_value(opt.final_value().as<bool>());
    } else if (opt.type() == typeid(double)) {
        vs = po::value<double>()
             ->default_value(opt.final_value().as<double>());
    } else if (opt.type() == typeid(std::string)) {
        po::typed_value<std::string>* tv = po::value<std::string>();
        tv->default_value(opt.final_value().as<std::string>());
        if (opt.required_) {
            tv->required();
        }
        vs = tv;
    } else if (opt.type() == typeid(std::vector<std::string>)) {
        po::typed_value<std::vector<std::string> >* tv;
        tv = po::value<std::vector<std::string> >()->multitoken();
        typedef std::vector<std::string> List;
        std::vector<std::string> list = opt.final_value().as<List>();
        using namespace boost::algorithm;
        std::string list_str = join(list, " ");
        tv->default_value(list, list_str);
        if (opt.required_) {
            tv->required();
        }
        vs = tv;
    }
    BOOST_ASSERT(vs);
    add_unique_options(desc)(name.c_str(), vs, opt.description_.c_str());
}

void Processor::add_options(po::options_description& desc) const {
    if (!no_options()) {
        // add self options
        po::options_description self_options_1;
        typedef Impl::Name2Option::value_type Pair;
        BOOST_FOREACH (const Pair& name_and_opt, impl_->opts_) {
            const Option& opt = name_and_opt.second;
            if (!is_ignored(opt.name_)) {
                add_option(self_options_1, opt_prefixed(opt.name_), opt);
            }
        }
        po::options_description self_options_2;
        add_options_impl(self_options_2);
        po::options_description self_options_2a;
        copy_not_ignored(self_options_2, self_options_2a);
        po::options_description new_opts(name());
        add_new_options(self_options_1, new_opts, &desc);
        add_new_options(self_options_2a, new_opts, &desc);
        if (!new_opts.options().empty()) {
            desc.add(new_opts);
        }
        BOOST_FOREACH (const Processor* child, impl_->children_) {
            child->add_options(desc);
        }
    }
}

void Processor::apply_options(const po::variables_map& vm0) {
    po::variables_map vm = vm0;
    if (no_options()) {
        // remove all options except --timing and --workers
        BOOST_FOREACH (const po::variables_map::value_type& key_value, vm0) {
            const std::string& key = key_value.first;
            if (key != opt_prefixed("timing") &&
                    key != opt_prefixed("workers")) {
                vm.erase(key);
            }
        }
    }
    // remove ignored options (even --timing and --workers can be ignored)
    typedef boost::shared_ptr<po::option_description> OptPtr;
    BOOST_FOREACH (OptPtr ignored_opt, impl_->ignored_options_.options()) {
        vm.erase(ignored_opt->long_name());
    }
    typedef Impl::Name2Option::value_type Pair;
    BOOST_FOREACH (const Pair& name_and_opt, impl_->opts_) {
        const std::string& name = name_and_opt.first;
        std::string prefixed_name = opt_prefixed(name);
        if (vm.count(prefixed_name)) {
            set_opt_value(name, vm[prefixed_name].value());
        }
    }
    apply_options_impl(vm);
    BOOST_FOREACH (Processor* child, impl_->children_) {
        child->apply_options(vm);
    }
}

std::vector<std::string> Processor::options_errors() const {
    std::vector<std::string> result;
    typedef Impl::Name2Option::value_type Pair;
    BOOST_FOREACH (const Pair& name_and_opt, impl_->opts_) {
        const Option& opt = name_and_opt.second;
        if (opt.required_) {
            if (opt.type() == typeid(std::string)) {
                if (opt.final_value().as<std::string>().empty()) {
                    result.push_back("Required option " + opt.name_ +
                                     " is empty");
                }
            }
            if (opt.type() == typeid(std::vector<std::string>)) {
                if (opt.final_value().as<std::vector<std::string> >().empty()) {
                    result.push_back("Required option " + opt.name_ +
                                     " is empty");
                }
            }
        }
    }
    BOOST_FOREACH (const OptionsChecker& checker, impl_->checkers_) {
        std::string message;
        bool valid = checker(message);
        if (!valid) {
            result.push_back(message);
        }
    }
    return result;
}

std::vector<std::string> Processor::options_warnings() const {
    std::vector<std::string> result;
    BOOST_FOREACH (const OptionsChecker& checker, impl_->checkers_) {
        std::string message;
        bool valid = checker(message);
        if (valid && !message.empty()) {
            result.push_back(message);
        }
    }
    return result;
}

void Processor::apply_vector_options(const std::vector<std::string>& options) {
    StringToArgv args;
    BOOST_FOREACH (const std::string& opt, options) {
        args.add_argument(opt);
    }
    po::options_description desc;
    add_options(desc);
    po::variables_map vm;
    po::store(po::command_line_parser(args.argc(), args.argv()).options(desc)
              .allow_unregistered().run(), vm);
    // po::notify(vm); // to pass required options check
    apply_options(vm);
}

void Processor::apply_string_options(const std::string& options) {
    using namespace boost::algorithm;
    using boost::tokenizer;
    using boost::escaped_list_separator;
    typedef tokenizer<escaped_list_separator<char> > tok_t;
    tok_t tok(options, escaped_list_separator<char>('\\', ' ', '\"'));
    std::vector<std::string> opts;
    BOOST_FOREACH (std::string opt, tok) {
        trim_right(opt);
        opts.push_back(opt);
    }
    apply_vector_options(opts);
}

bool Processor::run() const {
    bool result = false;
    if (workers() != 0 && block_set()) {
        using namespace boost::posix_time;
        ptime before, after;
        if (timing()) {
            before = microsec_clock::universal_time();
        }
        result = run_impl();
        if (timing()) {
            after = microsec_clock::universal_time();
            impl_->milliseconds_ += (after - before).total_milliseconds();
            key(); // to memorize value. RTTI would be invalid in ~Processor()
        }
    }
    return result;
}

bool Processor::apply_to_block(Block* block) const {
    return apply_to_block_impl(block);
}

std::string Processor::name() const {
    if (!impl_->name_.empty()) {
        return impl_->name_.c_str();
    } else {
        const char* ni = name_impl();
        if (*ni == '\0') {
            // empty
            return key();
        } else {
            return ni;
        }
    }
}

void Processor::set_name(const std::string& name) {
    impl_->name_ = name;
}

bool Processor::apply(const BlockSetPtr& bs) const {
    BlockSetPtr prev = block_set();
    const_cast<Processor*>(this)->set_block_set(bs);
    bool result = run();
    const_cast<Processor*>(this)->set_block_set(prev);
    return result;
}

std::string Processor::key() const {
    if (impl_->key_.empty()) {
        impl_->key_ = processor_name(this);
    }
    return impl_->key_;
}

void Processor::set_key(const std::string& key) {
    impl_->key_ = key;
}

Processor* Processor::parent() const {
    return impl_->parent_;
}

static void remove_child(std::vector<Processor*>& children, Processor* child) {
    children.erase(std::remove(children.begin(), children.end(), child),
                   children.end()); // TODO template remove_from_vector()
}

void Processor::set_parent(Processor* parent) {
    if (parent != impl_->parent_) {
        if (impl_->parent_) {
            remove_child(impl_->parent_->impl_->children_, this);
        }
        impl_->parent_ = parent;
        if (parent) {
            parent->impl_->children_.push_back(this);
        }
    }
}

std::vector<Processor*> Processor::children() const {
    return impl_->children_;
}

Meta* Processor::meta() const {
    if (impl_->meta_) {
        return impl_->meta_;
    } else if (parent()) {
        return parent()->meta();
    } else {
        return tss_meta();
    }
}

void Processor::set_meta(Meta* meta) {
    impl_->meta_ = meta;
}

const std::string& Processor::opt_prefix() const {
    return impl_->opt_prefix_;
}

void Processor::set_opt_prefix(const std::string& opt_prefix) {
    impl_->opt_prefix_ = opt_prefix;
}

std::string Processor::opt_prefixed(const std::string& name) const {
    std::vector<const Processor*> ancestors;
    const Processor* processor = this;
    while (processor) {
        ancestors.push_back(processor);
        processor = processor->parent();
    }
    std::string result;
    BOOST_REVERSE_FOREACH (const Processor* p, ancestors) {
        result += p->opt_prefix();
    }
    result += name;
    return result;
}

std::vector<std::string> Processor::opts() const {
    std::vector<std::string> result;
    typedef Impl::Name2Option::value_type Pair;
    BOOST_FOREACH (const Pair& name_and_opt, impl_->opts_) {
        const Option& opt = name_and_opt.second;
        result.push_back(opt.name_);
    }
    return result;
}

bool Processor::has_opt(const std::string& name) const {
    return impl_->opts_.find(name) != impl_->opts_.end();
}

const std::string& Processor::opt_description(const std::string& name) const {
    typedef Impl::Name2Option::const_iterator It;
    It it = impl_->opts_.find(name);
    if (it == impl_->opts_.end()) {
        throw Exception("No option with name '" + name + "'");
    } else {
        return it->second.description_;
    }
}

const std::type_info& Processor::opt_type(const std::string& name) const {
    typedef Impl::Name2Option::const_iterator It;
    It it = impl_->opts_.find(name);
    if (it == impl_->opts_.end()) {
        throw Exception("No option with name '" + name + "'");
    } else {
        return it->second.type();
    }
}

const AnyAs& Processor::default_opt_value(const std::string& name) const {
    typedef Impl::Name2Option::const_iterator It;
    It it = impl_->opts_.find(name);
    if (it == impl_->opts_.end()) {
        throw Exception("No option with name '" + name + "'");
    } else {
        return it->second.default_value_;
    }
}

const AnyAs& Processor::opt_value(const std::string& name) const {
    typedef Impl::Name2Option::const_iterator It;
    It it = impl_->opts_.find(name);
    if (it == impl_->opts_.end()) {
        throw Exception("No option with name '" + name + "'");
    } else {
        return it->second.final_value();
    }
}

void Processor::set_opt_value(const std::string& name,
                              const AnyAs& value) {
    typedef Impl::Name2Option::iterator It;
    It it = impl_->opts_.find(name);
    if (it == impl_->opts_.end()) {
        throw Exception("No option with name '" + name + "'");
    }
    Option& opt = it->second;
    AnyAs v = value;
    BOOST_FOREACH (const OptionValidator& validator, opt.validators_) {
        v = validator(v);
    }
    if (!v.empty() && opt.type() != v.type()) {
        throw Exception("Type of value of option "
                        "'" + name + "' (" + v.type().name() + ") "
                        "differs from type of default value "
                        "(" + opt.type().name() + ")");
    }
    opt.value_ = v;
}

void Processor::add_options_impl(po::options_description& desc) const
{ }

void Processor::apply_options_impl(const po::variables_map& vm)
{ }

bool Processor::run_impl() const {
    return false;
}

bool Processor::apply_to_block_impl(Block* block) const {
    BlockSetPtr block_set = new_bs();
    block_set->insert(block);
    bool result = apply(block_set);
    block_set->detach(block);
    return result;
}

const char* Processor::name_impl() const {
    return "";
}

void Processor::add_opt(const std::string& name,
                        const std::string& description,
                        const AnyAs& default_value,
                        bool required) {
    BOOST_ASSERT(good_opt_type(default_value.type()));
    impl_->opts_[name] = Option(name, description, default_value, required);
}

void Processor::remove_opt(const std::string& name, bool apply_prefix) {
    impl_->opts_.erase(apply_prefix ? opt_prefixed(name) : name);
}

static bool any_equal(const AnyAs& a, const AnyAs& b) {
    BOOST_ASSERT(good_opt_type(a.type()));
    BOOST_ASSERT(good_opt_type(b.type()));
    if (a.type() != b.type()) {
        return false;
    }
    if (a.type() == typeid(int)) {
        return a.as<int>() == b.as<int>();
    } else if (a.type() == typeid(bool)) {
        return a.as<bool>() == b.as<bool>();
    } else if (a.type() == typeid(double)) {
        return a.as<double>() == b.as<double>();
    } else if (a.type() == typeid(std::string)) {
        return a.as<std::string>() == b.as<std::string>();
    } else if (a.type() == typeid(std::vector<std::string>)) {
        return a.as<std::vector<std::string> >() ==
               b.as<std::vector<std::string> >();
    }
    throw Exception("wrong type of any");
}

void Processor::add_opt_validator(const std::string& name,
                                  const OptionValidator& validator) {
    BOOST_ASSERT(has_opt(name));
    Option& opt = impl_->opts_[name];
    BOOST_ASSERT(any_equal(validator(opt.default_value_), opt.default_value_));
    impl_->opts_[name].validators_.push_back(validator);
}

void Processor::add_opt_check(const OptionsChecker& checker) {
    impl_->checkers_.push_back(checker);
}

static double double_option(const Processor* p, const std::string& name) {
    if (p->opt_type(name) == typeid(double)) {
        return p->opt_value(name).as<double>();
    } else if (p->opt_type(name) == typeid(int)) {
        return double(p->opt_value(name).as<int>());
    } else {
        throw Exception("Bad option type (" +
                        std::string(p->opt_type(name).name()) +
                        "), must be int or double");
    }
}

static double double_transparent(double value) {
    return value;
}

static bool general_checker(bool result, const std::string& s, std::string& d) {
    if (!result) {
        d = s;
    }
    return result;
}

void Processor::add_opt_rule(const std::string& rule,
                             const std::string& message) {
    using namespace boost::algorithm;
    std::vector<std::string> parts;
    split(parts, rule, isspace, token_compress_on);
    const std::string& left_opt = parts[0];
    const std::string& op = parts[1];
    const std::string& right = parts[2];
    if (!has_opt(left_opt)) {
        throw Exception("No such option: " + left_opt);
    }
    const std::type_info& left_opt_type = impl_->opts_[left_opt].type();
    if (left_opt_type != typeid(int) && left_opt_type != typeid(double)) {
        throw Exception("Option type for rule must be int or double, not " +
                        std::string(left_opt_type.name()));
    }
    typedef boost::function<double()> Getter;
    Getter left_getter = boost::bind(double_option, this, left_opt);
    boost::function<bool(double, double)> cmp;
    if (op == "<") {
        cmp = std::less<double>();
    } else if (op == ">") {
        cmp = std::greater<double>();
    } else if (op == "<=") {
        cmp = std::less_equal<double>();
    } else if (op == ">=") {
        cmp = std::greater_equal<double>();
    } else {
        throw Exception("Operators for rule must be <.>,<=,>=, not " + op);
    }
    Getter right_getter;
    if (has_opt(right)) {
        const std::type_info& right_opt_type = impl_->opts_[right].type();
        if (right_opt_type != typeid(int) && right_opt_type != typeid(double)) {
            throw Exception("Option type for rule must be int or double, not " +
                            std::string(right_opt_type.name()));
        }
        right_getter = boost::bind(double_option, this, right);
    } else {
        // may throw here if bad value
        double right_value = boost::lexical_cast<double>(right);
        right_getter = boost::bind(double_transparent, right_value);
    }
    OptionsChecker checker = boost::bind(general_checker, boost::bind(cmp,
                                         boost::bind(left_getter),
                                         boost::bind(right_getter)),
                                         message, _1);
    add_opt_check(checker);
}

void Processor::add_opt_rule(const std::string& rule) {
    add_opt_rule(rule, rule);
}

void Processor::log_processor(std::ostream& o, int depth) {
    using namespace boost::posix_time;
    impl_->logged_ = true;
    const int TAB_SIZE = 4;
    o << std::string(depth * TAB_SIZE, ' '); // indent
    o << key() + ": ";
    o << to_simple_string(milliseconds(impl_->milliseconds_));
    o << std::endl;
    BOOST_FOREACH (Processor* child, impl_->children_) {
        child->log_processor(o, depth + 1);
    }
}

void Processor::copy_not_ignored(const po::options_description& source,
                                 po::options_description& dest) const {
    typedef boost::shared_ptr<po::option_description> OptPtr;
    BOOST_FOREACH (OptPtr opt, source.options()) {
        bool good_option = true;
        const Processor* p = this;
        while (p) {
            if (p->impl_->ignored_options_.find_nothrow(opt->long_name(),
                    /* approx */ false)) {
                good_option = false;
                break;
            }
            p = p->parent();
        }
        if (good_option) {
            dest.add(opt);
        }
    }
}

std::string processor_name(const Processor* processor) {
    return class_name(typeid(*processor).name());
}

}

