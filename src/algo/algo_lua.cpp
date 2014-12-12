/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <memory>
#include <set>
#include <sstream>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "luabind-error.hpp"
#include "luabind-format-signature.hpp"
#include <luabind/luabind.hpp>
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
#include "AlignmentRow.hpp"
#include "FragmentCollection.hpp"
#include "Pipe.hpp"
#include "BlocksJobs.hpp"
#include "Meta.hpp"
#include "AbstractAligner.hpp"
#include "cast.hpp"

namespace luabind {

template <typename T>
object type_to_lua(lua_State* L, const T& v) {
    default_converter<T>().to(L, v);
    object o = object(from_stack(L, -1));
    lua_pop(L, 1);
    return o;
}

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

typedef default_converter<npge::MapAny> dcM;

int dcM::compute_score(lua_State* L, int index) {
    if (lua_type(L, index) != LUA_TTABLE) {
        return -1;
    }
    // make sure that type of key is string
    lua_pushnil(L);
    if (index < 0) {
        index -= 1;
    }
    if (lua_next(L, index) != 0) {
        luabind::object o(from_stack(L, -1));
        int type = lua_type(L, -2);
        lua_pop(L, 2);
        return (type == LUA_TSTRING) ? 0 : -1;
    } else {
        // empty map
        return 0;
    }
}

template<typename F>
void for_all_types(F f) {
    f.template apply_type<npge::EmptyTable>();
    f.template apply_type<bool>();
    f.template apply_type<int>();
    f.template apply_class<npge::Decimal>();
    f.template apply_type<std::string>();
    f.template apply_type<npge::Strings>();
    f.template apply_class<boost::shared_ptr<npge::SetFc> >();
    f.template apply_class<boost::shared_ptr<npge::VectorFc> >();
    f.template apply_class<npge::SequencePtr>();
    f.template apply_class<npge::Fragment*>();
    f.template apply_class<npge::AlignmentRow*>();
    f.template apply_class<npge::Block*>();
    f.template apply_class<npge::BlockSetPtr>();
    f.template apply_class<npge::Processor*>();
    f.template apply_type<npge::Blocks>();
    f.template apply_type<npge::Fragments>();
}

struct AnyFrom {
    npge::AnyAs& a_;
    const object& o_;

    AnyFrom(npge::AnyAs& a, const object& o):
        a_(a), o_(o) {
    }

    template<typename T>
    void apply() const {
        if (a_.empty()) {
            try {
                a_ = object_cast<T>(o_);
            } catch (...) {
            }
        }
    }

    template<typename T>
    void apply_type() const {
        apply<T>();
    }

    template<typename T>
    void apply_class() const {
        apply<T>();
    }
};

npge::MapAny dcM::from(lua_State* L, int index) {
    using namespace npge;
    object table(from_stack(L, index));
    ASSERT_EQ(luabind::type(table), LUA_TTABLE);
    MapAny result;
    for (luabind::iterator it(table), end; it != end; ++it) {
        std::string key = object_cast<std::string>(it.key());
        AnyAs value;
        for_all_types(AnyFrom(value, *it));
        if (!value.empty()) {
            result[key] = value;
        }
    }
    return result;
}

struct AnyTo {
    const npge::AnyAs& a_;
    object& o_;
    lua_State* L_;

    AnyTo(const npge::AnyAs& a, object& o, lua_State* L):
        a_(a), o_(o), L_(L) {
    }

    template<typename T>
    void apply_class() const {
        if (!o_) {
            try {
                o_ = object(L_, a_.as<T>());
            } catch (...) {
            }
        }
    }

    template<typename T>
    void apply_type() const {
        if (!o_) {
            try {
                o_ = type_to_lua<T>(L_, a_.as<T>());
            } catch (...) {
            }
        }
    }
};

void dcM::to(lua_State* L, const npge::MapAny& a) {
    using namespace npge;
    lua_createtable(L, 0, a.size());
    BOOST_FOREACH (const MapAny::value_type& kv, a) {
        const std::string& key = kv.first;
        const AnyAs& value = kv.second;
        luabind::object o;
        for_all_types(AnyTo(value, o, L));
        if (o) {
            lua_pushstring(L, key.c_str());
            o.push(L);
            lua_rawset(L, -3);
        }
    }
}

}

namespace npge {

static std::string dumpf(const luabind::object& func) {
    using namespace luabind;
    object dump = globals(func.interpreter())["string"]["dump"];
    return object_cast<std::string>(dump(func));
}

static luabind::object loads(lua_State* L,
                             const std::string& f) {
    using namespace luabind;
    object ls = globals(L)["loadstring"];
    return ls(f);
}

static Processor* new_processor() {
    return new Processor;
}

static void delete_processor(Processor* p) {
    delete p;
}

static SharedProcessor get_processor_deleter(Processor* p) {
    return SharedProcessor(p);
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
    print_help(":stdout", p);
}

static void processor_print_tree(Processor* p) {
    print_processor_tree(":stdout", p);
}

static luabind::scope register_processor() {
    using namespace luabind;
    return class_<Processor>("Processor")
           .scope [
               def("new", &new_processor),
               def("delete", &delete_processor),
               def("deleter", &get_processor_deleter)
           ]
           .def(tostring(self))
           .def("get_block_sets", &processor_get_block_sets)
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
           .def("print_tree", &processor_print_tree)
          ;
}

class LuaProcessor: public Processor {
public:
    void set_action(const luabind::object& f) {
        f_ = dumpf(f);
    }

protected:
    void run_impl() const {
        if (!f_.empty()) {
            luabind::object f = loads(meta()->L(), f_);
            f(const_cast<LuaProcessor*>(this));
        }
    }

private:
    std::string f_;
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
           .def("max_iterations", &Pipe::max_iterations)
           .def("set_max_iterations", &Pipe::set_max_iterations)
           .def("processors", &Pipe::processors)
          ;
}

struct LuaWD : public WorkData {
    luabind::object work_o_;
    MapAny work_a_;
};

struct LuaTD : public ThreadData {
    luabind::object work_o_;
    luabind::object thread_o_;
    MapAny thread_a_;
};

class LuaBlocksJobs : public BlocksJobs {
public:
    void set_change_blocks(const luabind::object& f) {
        change_blocks_ = f;
    }

    void change_blocks_impl(Blocks& bb) const {
        using namespace luabind;
        if (change_blocks_) {
            bb = object_cast<Blocks>(change_blocks_(bb));
        }
    }

    void set_before_work(const luabind::object& f) {
        before_work_ = f;
    }

    LuaWD* before_work_impl() const {
        using namespace luabind;
        std::auto_ptr<LuaWD> wd(new LuaWD);
        if (before_work_) {
            wd->work_o_ = before_work_();
        } else {
            wd->work_o_ = luabind::newtable(meta()->L());
        }
        wd->work_o_["processor"] = (Processor*)this;
        try {
            wd->work_a_ = object_cast<MapAny>(wd->work_o_);
        } catch (...) {
        }
        return wd.release();
    }

    ThreadData* before_thread_impl() const {
        return new LuaTD;
    }

    void set_initialize_thread(const luabind::object& f) {
        initialize_thread_ = dumpf(f);
    }

    void initialize_thread_impl(ThreadData* d0) const {
        using namespace luabind;
        ASSERT_TRUE(d0);
        LuaTD* td = D_CAST<LuaTD*>(d0);
        ASSERT_TRUE(td);
        LuaWD* wd = D_CAST<LuaWD*>(td->work_data());
        ASSERT_TRUE(wd);
        td->work_o_ = type_to_lua<MapAny>(meta()->L(),
                                          wd->work_a_);
        if (!initialize_thread_.empty()) {
            object f = loads(meta()->L(), initialize_thread_);
            td->thread_o_ = f(td->work_o_);
        } else {
            td->thread_o_ = luabind::newtable(meta()->L());
        }
    }

    void set_process_block(const luabind::object& f) {
        process_block_ = dumpf(f);
    }

    void process_block_impl(Block* block,
                            ThreadData* d0) const {
        using namespace luabind;
        if (!process_block_.empty()) {
            object f = loads(meta()->L(), process_block_);
            LuaTD* td = D_CAST<LuaTD*>(d0);
            ASSERT_TRUE(td);
            f(block, td->thread_o_, td->work_o_);
        }
    }

    void finish_thread_impl(ThreadData* d0) const {
        using namespace luabind;
        ASSERT_TRUE(d0);
        LuaTD* td = D_CAST<LuaTD*>(d0);
        ASSERT_TRUE(td);
        try {
            td->thread_a_ = object_cast<MapAny>(td->thread_o_);
        } catch (...) {
        }
    }

    void set_after_thread(const luabind::object& f) {
        after_thread_ = f;
    }

    void after_thread_impl(ThreadData* d0) const {
        using namespace luabind;
        if (after_thread_) {
            ASSERT_TRUE(d0);
            LuaTD* td = D_CAST<LuaTD*>(d0);
            ASSERT_TRUE(td);
            object thread_o = type_to_lua<MapAny>(meta()->L(),
                                                  td->thread_a_);
            LuaWD* wd = D_CAST<LuaWD*>(td->work_data());
            ASSERT_TRUE(wd);
            after_thread_(thread_o, wd->work_o_);
        }
    }

    void set_after_work(const luabind::object& f) {
        after_work_ = f;
    }

    void after_work_impl(WorkData* wd0) const {
        using namespace luabind;
        if (after_work_) {
            LuaWD* wd = D_CAST<LuaWD*>(wd0);
            ASSERT_TRUE(wd);
            after_work_(wd->work_o_);
        }
    }

private:
    mutable luabind::object change_blocks_;
    mutable luabind::object before_work_;
    mutable std::string initialize_thread_;
    mutable std::string process_block_;
    mutable luabind::object after_thread_;
    mutable luabind::object after_work_;
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
           .def("set_before_work",
                &LuaBlocksJobs::set_before_work)
           .def("set_initialize_thread",
                &LuaBlocksJobs::set_initialize_thread)
           .def("set_process_block",
                &LuaBlocksJobs::set_process_block)
           .def("set_after_thread",
                &LuaBlocksJobs::set_after_thread)
           .def("set_after_work",
                &LuaBlocksJobs::set_after_work)
          ;
}

static Processor* return_processor(Meta* meta, std::string f,
                                   std::string key) {
    using namespace luabind;
    object func = loads(meta->L(), f);
    Processor* p = object_cast<Processor*>(func());
    if (p) {
        p->set_key(key);
    }
    return p;
}

static void meta_set_returner2(
    Meta* meta,
    const luabind::object& func,
    const std::string& key) {
    std::string f = dumpf(func);
    bool overwrite = true;
    meta->set_returner(boost::bind(return_processor,
                                   meta, f, key),
                       key, overwrite);
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
    print_config(":stdout", meta);
}

static luabind::scope register_meta() {
    using namespace luabind;
    return class_<Meta>("Meta")
           .def("has", &Meta::has)
           .def("get_plain", &Meta::get_plain)
           .def("set_returner", &meta_set_returner2)
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
           .def("get_section", &Meta::get_section)
           .def("set_section", &Meta::set_section)
           .def("set_opt", &Meta::set_opt)
           .def("set_opt", &meta_set_opt)
           .def("set_opt", &meta_set_opt1)
           .def("set_opt", &meta_set_opt2)
           .def("set_opt_func", &meta_set_opt_func)
           .def("has_opt", &Meta::has_opt)
           .def("opts", &Meta::opts)
           .def("sections", &Meta::sections)
           .def("opts_of_section", &Meta::opts_of_section)
           .def("remove_opt", &Meta::remove_opt)
           .def("print_config", &meta_print_config)
           .def("print_config", &meta_print_config1)
          ;
}

}

extern "C" int init_algo_lua(lua_State* L) {
    using namespace luabind;
    using namespace npge;
    luabind::open(L);
    module(L) [
        register_processor(),
        register_luaprocessor(),
        register_pipe(),
        register_blocks_jobs(),
        register_meta(),
        def("remove_pure_gap_columns",
            &AbstractAligner::remove_pure_gap_columns)
    ];
    return 0;
}

