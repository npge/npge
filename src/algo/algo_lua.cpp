/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <boost/bind.hpp>

#include <luabind/luabind.hpp>
#include <luabind/operator.hpp>
#include <luabind/object.hpp>

#include "algo_lua.hpp"
#include "model_lua.hpp"
#include "util_lua.hpp"
#include "global.hpp"
#include "Processor.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"

namespace npge {

static Processor* new_processor() {
    return new Processor;
}

static void delete_processor(Processor* p) {
    delete p;
}

static void processor_set_options(Processor* p,
                                  const std::string& options) {
    p->set_options(options);
}

static Strings processor_get_block_sets(Processor* p) {
    Strings result;
    p->get_block_sets(result);
    return result;
}

static void processor_add_opt_rule0(Processor* p,
                                    const std::string& rule) {
    p->add_opt_rule(rule);
}

static void processor_add_opt_rule1(Processor* p,
                                    const std::string& rule,
                                    const std::string& m) {
    p->add_opt_rule(rule, m);
}

static AnyAs processor_go(Processor* p,
                          const std::string& key) {
    return p->go(key);
}

static void processor_add_opt(Processor* p,
                              const std::string& name,
                              const std::string& descr,
                              const AnyAs& value) {
    p->add_opt(name, descr, value);
}

static void processor_add_gopt(Processor* p,
                               const std::string& name,
                               const std::string& descr,
                               const std::string& gon) {
    p->add_gopt(name, descr, gon);
}

static void processor_remove_opt(Processor* p,
                                 const std::string& name) {
    p->remove_opt(name);
}

static AnyAs validate(luabind::object f,
                      const AnyAs& v) {
    using namespace npge;
    using namespace luabind;
    if (v.type() == typeid(bool)) {
        return object_cast<bool>(f(v.as<bool>()));
    } else if (v.type() == typeid(int)) {
        return object_cast<int>(f(v.as<int>()));
    } else if (v.type() == typeid(Decimal)) {
        // TODO
        return v;
    } else if (v.type() == typeid(std::string)) {
        return object_cast<std::string>(f(v.as<std::string>()));
    } else if (v.type() == typeid(Strings)) {
        return object_cast<Strings>(f(v.as<Strings>()));
    }
    // wrong type
    return v;
}

static void processor_add_opt_validator(
    Processor* p,
    const std::string& name,
    const luabind::object& validator) {
    p->add_opt_validator(name, boost::bind(validate,
                                           validator, _1));
}

static std::string processor_opt_type(Processor* p,
                                      const std::string& name) {
    return p->default_opt_value(name).type_name();
}

struct ProcessorWrapper : public Processor, luabind::wrap_base {
    virtual void run_impl() {
        call<void>("run_impl");
    }

    void public_run_impl() {
        run_impl();
    }

    void default_run_impl() {
        Processor::run_impl();
    }
};

static luabind::scope register_processor() {
    using namespace luabind;
    return class_<Processor, ProcessorWrapper>("Processor")
           .scope [
               def("new", &new_processor),
               def("delete", &delete_processor)
           ]
           .def(tostring(self))
           .def("declare_bs", &Processor::declare_bs)
           .def("remove_bs", &Processor::remove_bs)
           .def("bs_description", &Processor::bs_description)
           .def("get_bs", &Processor::get_bs)
           .def("set_bs", &Processor::set_bs)
           .def("has_bs", &Processor::has_bs)
           .def("point_bs", &Processor::point_bs)
           .def("set_options", &Processor::set_options)
           .def("set_options", &processor_set_options)
           .def("block_set", &Processor::block_set)
           .def("set_block_set", &Processor::set_block_set)
           .def("other", &Processor::other)
           .def("set_other", &Processor::set_other)
           .def("get_block_sets", &processor_get_block_sets)
           .def("workers", &Processor::workers)
           .def("set_workers", &Processor::set_workers)
           .def("write_log", &Processor::write_log)
           .def("close_log", &Processor::close_log)
           .def("no_options", &Processor::no_options)
           .def("set_no_options", &Processor::set_no_options)
           .def("add_ignored_option",
                &Processor::add_ignored_option)
           .def("is_ignored", &Processor::is_ignored)
           .def("timing", &Processor::timing)
           .def("set_timing", &Processor::set_timing)
           .def("assign", &Processor::assign)
           // TODO add_opt_check
           .def("add_opt_rule", &processor_add_opt_rule0)
           .def("add_opt_rule", &processor_add_opt_rule1)
           .def("options_errors", &Processor::options_errors)
           .def("options_warnings",
                &Processor::options_warnings)
           .def("apply_vector_options",
                &Processor::apply_vector_options)
           .def("apply_string_options",
                &Processor::apply_string_options)
           .def("run", &Processor::run)
           .def("apply", &Processor::apply)
           .def("apply_to_block", &Processor::apply_to_block)
           .def("name", &Processor::name)
           .def("set_name", &Processor::set_name)
           .def("key", &Processor::key)
           .def("set_key", &Processor::set_key)
           .def("parent", &Processor::parent)
           .def("set_parent", &Processor::set_parent)
           // TODO children
           .def("clone", &Processor::clone)
           // TODO meta, set_meta
           .def("go", &Processor::go)
           .def("go", &processor_go)
           .def("opt_prefix", &Processor::opt_prefix)
           .def("set_opt_prefix", &Processor::set_opt_prefix)
           .def("opt_prefixed", &Processor::opt_prefixed)
           .def("add_opt", &Processor::add_opt)
           .def("add_opt", &processor_add_opt)
           .def("add_gopt", &Processor::add_gopt)
           .def("add_gopt", &processor_add_gopt)
           .def("remove_opt", &Processor::remove_opt)
           .def("remove_opt", &processor_remove_opt)
           .def("add_opt_validator",
                &processor_add_opt_validator)
           .def("opts", &Processor::opts)
           .def("has_opt", &Processor::has_opt)
           .def("opt_description", &Processor::opt_description)
           .def("opt_type", &processor_opt_type)
           .def("default_opt_value",
                &Processor::default_opt_value)
           .def("opt_value", &Processor::opt_value)
           .def("set_opt_value", &Processor::set_opt_value)
           // TODO set_opt_getter
           .def("fix_opt_value", &Processor::fix_opt_value)
           // TODO fix_opt_getter
           .def("interrupt", &Processor::interrupt)
           .def("is_interrupted", &Processor::is_interrupted)
           .def("tmp_file", &Processor::tmp_file)
           .def("processor_name", &processor_name)
           .def("run_impl", &ProcessorWrapper::public_run_impl,
                &ProcessorWrapper::default_run_impl)
          ;
}

}

extern "C" int init_algo_lua(lua_State* L) {
    using namespace luabind;
    using namespace npge;
    open(L);
    module(L) [
        register_processor()
    ];
    return 0;
}

