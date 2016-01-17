/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>

#include "Processor.hpp"
#include "Filter.hpp"
#include "Pipe.hpp"
#include "Filter.hpp"
#include "Decimal.hpp"

using namespace npge;

BOOST_AUTO_TEST_CASE (processor_name_main) {
    SharedProcessor p(new Filter);
    BOOST_CHECK(processor_name(p.get()) == "Filter");
    SharedProcessor p1(new Filter);
    BOOST_CHECK(processor_name(p1.get()) == "Filter");
}

BOOST_AUTO_TEST_CASE (processor_set_options) {
    Processor p1;
    p1.set_timing(false);
    BOOST_CHECK(!p1.no_options());
    BOOST_CHECK(!p1.timing());
    //
    Processor p2;
    p2.set_options("no_options");
    BOOST_CHECK(p2.no_options());
    //
    Processor p3;
    p3.set_options("--timing=1");
    BOOST_CHECK(p3.timing());
    p3.set_timing(false); // not to write to std err
}

class NoOptionsPipe : public Pipe {
public:
    Filter* or2_;

    NoOptionsPipe() {
        or2_ = new Filter;
        or2_->set_opt_value("min-fragment", 10);
        add(or2_);
    }
};

BOOST_AUTO_TEST_CASE (processor_NoOptionsPipe) {
    NoOptionsPipe nop;
    nop.set_options("--min-fragment=10");
    BOOST_CHECK(nop.or2_->opt_value("min-fragment").as<int>() == 10);
    nop.apply_string_options("--min-fragment=20");
    BOOST_CHECK(nop.or2_->opt_value("min-fragment").as<int>() == 20);
    nop.set_no_options(true);
    nop.apply_string_options("--min-fragment=30");
    BOOST_CHECK(nop.or2_->opt_value("min-fragment").as<int>() == 20);
    nop.set_options("--min-fragment=10");
    BOOST_CHECK(nop.or2_->opt_value("min-fragment").as<int>() == 20);
}

#define S std::string

class TestProcessor : public Processor {
public:
    TestProcessor() {
        add_opt("bool", "Boolean", false);
        add_opt("integer", "Integer for testing", 42);
        add_opt("decimal", "Decimal", D(1.46));
        add_opt("string", "string", S(""), /* required */ false);
        add_opt("string-0", "Null string", S(""), /* required */ true);
        add_opt("string-1", "First string", S("first"), /* required */ true);
        add_opt("string-2", "Second string", S("default value of second"));
        Strings list;
        list.push_back("1");
        list.push_back("2");
        add_opt("list", "List of strings", list);
        add_opt("list-0", "", Strings(), false);
        add_opt("list-1", "", Strings(), true);
    }
};

BOOST_AUTO_TEST_CASE (processor_options_main) {
    TestProcessor p;
    BOOST_CHECK(!p.opts().empty());
    BOOST_FOREACH (std::string opt_name, p.opts()) {
        BOOST_CHECK(p.has_opt(opt_name));
    }
    BOOST_CHECK(p.has_opt("bool"));
    BOOST_CHECK(p.opt_type("bool") == typeid(bool));
    BOOST_CHECK(p.opt_description("bool") == "Boolean");
    BOOST_CHECK(!p.has_opt("bar"));
    BOOST_CHECK(p.opt_value("bool").as<bool>() == false);
    BOOST_CHECK(p.default_opt_value("bool").as<bool>() == false);
    p.set_opt_value("bool", true);
    BOOST_CHECK(p.opt_value("bool").as<bool>() == true);
    BOOST_CHECK(p.default_opt_value("bool").as<bool>() == false);
    BOOST_CHECK(p.opt_value("integer").as<int>() == 42);
    bool error_caught = false;
    try {
        p.opt_value("integer").as<bool>();
    } catch (...) {
        error_caught = true;
    }
    BOOST_CHECK(error_caught);
    BOOST_CHECK(p.opt_value("string").as<std::string>() == "");
    p.set_opt_value("string", S("test"));
    BOOST_CHECK(p.opt_value("string").as<std::string>() == "test");
}

BOOST_AUTO_TEST_CASE (processor_options_apply) {
    TestProcessor p;
    BOOST_CHECK(!p.options_errors().empty());
    p.apply_string_options("");
    BOOST_CHECK(!p.options_errors().empty());
    p.apply_string_options("--string-0='111' --list-1 l i s t");
    BOOST_CHECK(p.options_errors().empty());
}

BOOST_AUTO_TEST_CASE (processor_options_apply_prefix) {
    TestProcessor p;
    BOOST_CHECK(p.opt_prefix() == "");
    BOOST_CHECK(p.has_opt("bool"));
    BOOST_CHECK(p.opt_prefixed("bool") == "bool");
    p.set_opt_prefix("foo-");
    BOOST_CHECK(p.opt_prefix() == "foo-");
    BOOST_CHECK(p.has_opt("bool"));
    BOOST_CHECK(!p.has_opt("foo-bool"));
    BOOST_CHECK(p.opt_prefixed("bool") == "foo-bool");
    BOOST_CHECK(!p.options_errors().empty());
    p.apply_string_options("");
    BOOST_CHECK(!p.options_errors().empty());
    p.apply_string_options("--foo-string-0='111' --foo-list-1 l i s t");
    BOOST_CHECK(p.options_errors().empty());
}

class RulesProcessor : public Processor {
public:
    RulesProcessor() {
        add_opt("integer", "Integer for testing", 42);
        add_opt_rule("integer < 50", "Integer must be less than 50");
        add_opt_rule("integer >= 10");
        add_opt("decimal", "Decimal", D(1.46));
        add_opt_rule("decimal < integer");
        add_opt_rule("integer > decimal");
        add_opt_rule("decimal > -1.2");
        add_opt_rule("decimal < 20");
    }
};

BOOST_AUTO_TEST_CASE (processor_options_rules) {
    RulesProcessor p;
    BOOST_CHECK(p.options_errors().empty());
    p.apply_string_options("");
    BOOST_CHECK(p.options_errors().empty());
    p.set_opt_value("integer", 49);
    BOOST_CHECK(p.options_errors().empty());
    p.set_opt_value("integer", 50);
    BOOST_CHECK(p.options_errors().size() == 1);
    p.set_opt_value("integer", 10);
    BOOST_CHECK(p.options_errors().empty());
    p.set_opt_value("decimal", D(19.0));
    BOOST_CHECK(p.options_errors().size() == 2);
    p.set_opt_value("decimal", D(20.0));
    BOOST_CHECK(p.options_errors().size() == 3);
    p.set_opt_value("decimal", D(-1.5));
    BOOST_CHECK(p.options_errors().size() == 1);
    p.set_opt_value("integer", 9);
    BOOST_CHECK(p.options_errors().size() == 2);
}

BOOST_AUTO_TEST_CASE (processor_workers) {
    Processor p;
    p.set_opt_value("workers", -1);
    BOOST_CHECK(p.opt_value("workers").as<int>() > 0);
}

BOOST_AUTO_TEST_CASE (processor_parent) {
    Processor* parent = new Processor;
    Processor* child = new Processor;
    child->set_parent(parent);
    delete child;
    delete parent;
}

