/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <set>
#include <sstream>
#include <boost/bind.hpp>

#include <luabind/luabind.hpp>
#include <luabind/tag_function.hpp>
#include <luabind/operator.hpp>
#include <luabind/object.hpp>

#include "algo_lua.hpp"
#include "model_lua.hpp"
#include "util_lua.hpp"
#include "global.hpp"
#include "throw_assert.hpp"
#include "Processor.hpp"
#include "process.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Pipe.hpp"
#include "BlocksJobs.hpp"
#include "Meta.hpp"
#include "FragmentCollection.hpp"

namespace luabind {

typedef default_converter<Processors> dcP;

int dcP::compute_score(lua_State* L, int index) {
    return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
}

Processors dcP::from(lua_State* L, int index) {
    // not implemented
    return Processors();
}

void dcP::to(lua_State* L, const Processors& a) {
    lua_createtable(L, a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        npge::Processor* p = a[i];
        luabind::object o(L, p);
        o.push(L);
        lua_rawseti(L, -2, i + 1);
    }
}

}

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

static bool check(luabind::object f, std::string& m) {
    // lua function should return {status, message}
    using namespace npge;
    using namespace luabind;
    object r;
    r = f();
    m = object_cast<std::string>(r[2]);
    return object_cast<bool>(r[1]);
}

static void processor_add_opt_check(
    Processor* p,
    const luabind::object& f) {
    p->add_opt_check(boost::bind(check, f, _1));
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

static void processor_add_opt1(Processor* p,
                               const std::string& name,
                               const std::string& descr,
                               const Decimal& value,
                               bool required) {
    p->add_opt(name, descr, value, required);
}

static void processor_add_opt2(Processor* p,
                               const std::string& name,
                               const std::string& descr,
                               const Decimal& value) {
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
        return object_cast<Decimal>(f(v.as<Decimal>()));
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

static void processor_set_opt_value(Processor* p,
                                    const std::string& name,
                                    const Decimal& value) {
    p->set_opt_value(name, value);
}

static AnyAs getter(AnyAs d, luabind::object f) {
    using namespace npge;
    using namespace luabind;
    if (d.type() == typeid(bool)) {
        return object_cast<bool>(f());
    } else if (d.type() == typeid(int)) {
        return object_cast<int>(f());
    } else if (d.type() == typeid(Decimal)) {
        return object_cast<Decimal>(f());
    } else if (d.type() == typeid(std::string)) {
        return object_cast<std::string>(f());
    } else if (d.type() == typeid(Strings)) {
        return object_cast<Strings>(f());
    }
    // wrong type
    return AnyAs();
}

static void processor_set_opt_getter(
    Processor* p,
    const std::string& name,
    const luabind::object& f) {
    AnyAs d = p->default_opt_value(name);
    p->set_opt_getter(name, boost::bind(getter, d, f));
}

static void processor_fix_opt_value(Processor* p,
                                    const std::string& name,
                                    const Decimal& value) {
    p->fix_opt_value(name, value);
}

static void processor_fix_opt_getter(
    Processor* p,
    const std::string& name,
    const luabind::object& f) {
    AnyAs d = p->default_opt_value(name);
    p->fix_opt_getter(name, boost::bind(getter, d, f));
}

static void processor_print_help(Processor* p) {
    print_help(":cerr", p);
}

static luabind::scope register_processor() {
    using namespace luabind;
    return class_<Processor>("Processor")
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
           .def("add_opt_check", &processor_add_opt_check)
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
           .def("children", &Processor::children)
           .def("clone", &Processor::clone)
           .def("meta", &Processor::meta)
           .def("set_meta", &Processor::set_meta)
           .def("go", &Processor::go)
           .def("go", &processor_go)
           .def("opt_prefix", &Processor::opt_prefix)
           .def("set_opt_prefix", &Processor::set_opt_prefix)
           .def("opt_prefixed", &Processor::opt_prefixed)
           .def("add_opt", &Processor::add_opt)
           .def("add_opt", &processor_add_opt)
           .def("add_opt", &processor_add_opt1)
           .def("add_opt", &processor_add_opt2)
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
           .def("set_opt_value", &processor_set_opt_value)
           .def("set_opt_getter", &processor_set_opt_getter)
           .def("fix_opt_value", &Processor::fix_opt_value)
           .def("fix_opt_value", &processor_fix_opt_value)
           .def("fix_opt_getter", &processor_fix_opt_getter)
           .def("interrupt", &Processor::interrupt)
           .def("is_interrupted", &Processor::is_interrupted)
           .def("tmp_file", &Processor::tmp_file)
           .def("processor_name", &processor_name)
           .def("print_help", &processor_print_help)
          ;
}

class LuaProcessor: public Processor {
public:
    void set_action(const luabind::object& f) {
        f_ = f;
    }

protected:
    void run_impl() const {
        if (f_) {
            f_();
        }
    }

private:
    mutable luabind::object f_;
};

static LuaProcessor* new_luaprocessor() {
    return new LuaProcessor;
}

static void delete_luaprocessor(LuaProcessor* p) {
    delete p;
}

static luabind::scope register_luaprocessor() {
    using namespace luabind;
    return class_<LuaProcessor, Processor>("LuaProcessor")
           .scope [
               def("new", &new_luaprocessor),
               def("delete", &delete_luaprocessor)
           ]
           .def("set_action", &LuaProcessor::set_action)
          ;
}

static Pipe* new_pipe() {
    return new Pipe;
}

static void delete_pipe(Pipe* p) {
    delete p;
}

static void pipe_add(Pipe* pipe, Processor* p) {
    pipe->add(p);
}

static void pipe_add1(Pipe* pipe, const std::string& key) {
    ASSERT_TRUE(pipe->meta());
    pipe->add(pipe->meta()->get_plain(key));
}

static void pipe_add2(Pipe* pipe, const std::string& key,
                      const std::string& options) {
    ASSERT_TRUE(pipe->meta());
    pipe->add(pipe->meta()->get_plain(key), options);
}

static luabind::scope register_pipe() {
    using namespace luabind;
    return class_<Pipe, Processor>("Pipe")
           .scope [
               def("new", &new_pipe),
               def("delete", &delete_pipe)
           ]
           .def("add", &Pipe::add)
           .def("add", &pipe_add)
           .def("add", &pipe_add1)
           .def("add", &pipe_add2)
           .def("max_loops", &Pipe::max_loops)
           .def("set_max_loops", &Pipe::set_max_loops)
           .def("processors", &Pipe::processors)
          ;
}

class LuaBlocksJobs : public BlocksJobs {
public:
    void run_impl() const {
        // Lua code shiould run in one thread
        const_cast<LuaBlocksJobs*>(this)->set_workers(1);
        BlocksJobs::run_impl();
    }

    void set_change_blocks(const luabind::object& f) {
        change_blocks_ = f;
    }

    void change_blocks_impl(Blocks& bb) const {
        using namespace luabind;
        if (change_blocks_) {
            bb = object_cast<Blocks>(change_blocks_(bb));
        }
    }

    void set_initialize_work(const luabind::object& f) {
        initialize_work_ = f;
    }

    void initialize_work_impl() const {
        if (initialize_work_) {
            initialize_work_();
        }
    }

    void set_process_block(const luabind::object& f) {
        process_block_ = f;
    }

    void process_block_impl(Block* block, ThreadData*) const {
        if (process_block_) {
            process_block_(block);
        }
    }

    void set_finish_work(const luabind::object& f) {
        finish_work_ = f;
    }

    void finish_work_impl() const {
        if (finish_work_) {
            finish_work_();
        }
    }

private:
    mutable luabind::object change_blocks_;
    mutable luabind::object initialize_work_;
    mutable luabind::object process_block_;
    mutable luabind::object finish_work_;
};

static LuaBlocksJobs* new_blocks_jobs() {
    return new LuaBlocksJobs;
}

static void delete_blocks_jobs(LuaBlocksJobs* p) {
    delete p;
}

static luabind::scope register_blocks_jobs() {
    using namespace luabind;
    return class_<LuaBlocksJobs, Processor>("BlocksJobs")
           .scope [
               def("new", &new_blocks_jobs),
               def("delete", &delete_blocks_jobs)
           ]
           .def("set_change_blocks",
                &LuaBlocksJobs::set_change_blocks)
           .def("set_initialize_work",
                &LuaBlocksJobs::set_initialize_work)
           .def("set_process_block",
                &LuaBlocksJobs::set_process_block)
           .def("set_finish_work",
                &LuaBlocksJobs::set_finish_work)
          ;
}

static Processor* return_processor(luabind::object f) {
    using namespace luabind;
    return object_cast<Processor*>(f());
}

static void meta_set_returner1(
    Meta* meta,
    const luabind::object& f,
    const std::string& key,
    bool overwrite) {
    meta->set_returner(boost::bind(return_processor, f),
                       key, overwrite);
}

static void meta_set_returner2(
    Meta* meta,
    const luabind::object& f,
    const std::string& key) {
    meta->set_returner(boost::bind(return_processor, f), key);
}

static void meta_set_returner3(
    Meta* meta,
    const luabind::object& f) {
    meta->set_returner(boost::bind(return_processor, f));
}

static AnyAs meta_get_opt(Meta* meta, const std::string& key) {
    return meta->get_opt(key);
}

static std::string meta_get_description(
    Meta* meta, const std::string& key) {
    return meta->get_description(key);
}

static void meta_set_opt(
    Meta* meta,
    const std::string& key,
    const AnyAs& value) {
    return meta->set_opt(key, value);
}

static void meta_set_opt1(
    Meta* meta,
    const std::string& key,
    const Decimal& value) {
    return meta->set_opt(key, value);
}

static void meta_set_opt2(
    Meta* meta,
    const std::string& key,
    const Decimal& value,
    const std::string& description) {
    return meta->set_opt(key, value, description);
}

static AnyAs opt_func(luabind::object f) {
    using namespace luabind;
    object r = f();
    if (type(r) == LUA_TBOOLEAN) {
        return object_cast<bool>(r);
    } else if (type(r) == LUA_TNUMBER) {
        return object_cast<int>(r);
    } else if (type(r) == LUA_TUSERDATA) {
        return object_cast<Decimal>(r);
    } else if (type(r) == LUA_TSTRING) {
        return object_cast<std::string>(r);
    } else if (type(r) == LUA_TTABLE) {
        return object_cast<Strings>(r);
    }
    return AnyAs();
}

static void meta_set_opt_func(
    Meta* meta,
    const std::string& key,
    const luabind::object& f) {
    meta->set_opt_func(key, boost::bind(opt_func, f));
}

static void meta_print_config(Meta* meta,
        const std::string& fname) {
    print_config(fname, meta);
}

static void meta_print_config1(Meta* meta) {
    print_config(":cerr", meta);
}

static luabind::scope register_meta() {
    using namespace luabind;
    return class_<Meta>("Meta")
           .def("has", &Meta::has)
           .def("get_plain", &Meta::get_plain)
           .def("set_returner", &meta_set_returner1)
           .def("set_returner", &meta_set_returner2)
           .def("set_returner", &meta_set_returner3)
           .def("keys", &Meta::keys)
           .def("empty", &Meta::empty)
           .def("clear", &Meta::clear)
           .def("placeholder_processor",
                &Meta::placeholder_processor)
           .def("reset_placeholder_processor",
                &Meta::reset_placeholder_processor)
           .def("get_opt", &Meta::get_opt)
           .def("get_opt", &meta_get_opt)
           .def("get_description", &Meta::get_description)
           .def("get_description", &meta_get_description)
           .def("set_description", &Meta::set_description)
           .def("set_opt", &Meta::set_opt)
           .def("set_opt", &meta_set_opt)
           .def("set_opt", &meta_set_opt1)
           .def("set_opt", &meta_set_opt2)
           .def("set_opt_func", &meta_set_opt_func)
           .def("has_opt", &Meta::has_opt)
           .def("opts", &Meta::opts)
           .def("remove_opt", &Meta::remove_opt)
           .def("print_config", &meta_print_config)
           .def("print_config", &meta_print_config1)
          ;
}

typedef std::set<Fragment*, FragmentCompare> FragmentsSet;
typedef FragmentCollection<Fragment*, FragmentsSet> SetFc;
typedef FragmentCollection<Fragment*, Fragments> VectorFc;

template<typename T>
struct find_overlap_fragments {
    Fragments operator()(T* fc, Fragment* f) const {
        Fragments result;
        fc->find_overlap_fragments(result, f);
        return result;
    }
};

template<typename T>
luabind::scope register_fragment_collection(const char* name) {
    using namespace luabind;
    return class_<T>(name)
           .def(constructor<>())
           .def("add_fragment", &T::add_fragment)
           .def("remove_fragment", &T::remove_fragment)
           .def("add_block", &T::add_block)
           .def("remove_block", &T::remove_block)
           .def("add_bs", &T::add_bs)
           .def("remove_bs", &T::remove_bs)
           .def("prepare", &T::prepare)
           .def("clear", &T::clear)
           .def("has_overlap", &T::has_overlap)
           .def("block_has_overlap", &T::block_has_overlap)
           .def("bs_has_overlap", &T::bs_has_overlap)
           .def("find_overlap_fragments",
                tag_function<Fragments(T*, Fragment*)>(
                    find_overlap_fragments<T>()))
           // TODO find_overlaps
          ;
}

}

extern "C" int init_algo_lua(lua_State* L) {
    using namespace luabind;
    using namespace npge;
    open(L);
    module(L) [
        register_processor(),
        register_luaprocessor(),
        register_pipe(),
        register_blocks_jobs(),
        register_meta(),
        register_fragment_collection<SetFc>("SetFc"),
        register_fragment_collection<VectorFc>("VectorFc")
    ];
    return 0;
}

