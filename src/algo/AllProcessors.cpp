/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "AllProcessors.hpp"
#include "Meta.hpp"
#include "Processor.hpp"
#include "process.hpp"
#include "name_to_stream.hpp"
#include "read_file.hpp"
#include "global.hpp"
#include "throw_assert.hpp"

namespace npge {

static void print_block_sets(std::ostream& out,
                             Processor* p) {
    Strings block_sets;
    p->get_block_sets(block_sets);
    out << "<ul>";
    BOOST_FOREACH (const std::string& bs_name, block_sets) {
        std::string descr = p->bs_description(bs_name);
        if (!descr.empty()) {
            out << "<li>";
            out << "<u>" << bs_name << "</u>" << ": " << descr;
            out << "</li>" << "\n";
        }
    }
    out << "</ul>";
}

typedef std::map<std::string, Strings> Key2Strings;

static void add_key(Key2Strings& k2s, Processor* p) {
    std::string key = p->key();
    set_sstream("help");
    print_help("help", p);
    std::string contents = read_file("help");
    remove_stream("help");
    using namespace boost::algorithm;
    Strings& lines = k2s[key];
    split(lines, contents, is_any_of("\n"));
}

static std::string replace_multi_spaces(std::string line) {
    int length = -1;
    while (line.length() != length) {
        length = line.length();
        using namespace boost::algorithm;
        replace_all(line, "  ", " ");
    }
    return line;
}

typedef std::set<std::string> StringSet;

static bool bad_line(const StringSet& cmn,
                     const std::string& line) {
    return cmn.find(replace_multi_spaces(line)) != cmn.end() ||
           line.find("-timing ") != std::string::npos ||
           line.find("-workers ") != std::string::npos;
}

static Strings remove_empty_sections(const Strings& lines) {
    Strings result;
    for (int i = 0; i < lines.size() - 1; i++) {
        const std::string& l = lines[i];
        const std::string& next = lines[i + 1];
        using namespace boost::algorithm;
        if (!l.empty() && ends_with(l, ":") && next.empty()) {
            i++;
        } else {
            result.push_back(l);
        }
    }
    return result;
}

static void remove_common(Key2Strings& k2s) {
    ASSERT_FALSE(k2s.empty());
    StringSet common;
    BOOST_FOREACH (std::string line, k2s.begin()->second) {
        line = replace_multi_spaces(line);
        if (!line.empty()) {
            common.insert(line);
        }
    }
    BOOST_FOREACH (const Key2Strings::value_type& k_s, k2s) {
        const Strings& lines = k_s.second;
        StringSet new_common;
        BOOST_FOREACH (std::string line, lines) {
            line = replace_multi_spaces(line);
            if (!line.empty()) {
                if (common.find(line) != common.end()) {
                    new_common.insert(line);
                }
            }
        }
        common.swap(new_common);
    }
    BOOST_FOREACH (Key2Strings::value_type& k_s, k2s) {
        Strings& lines = k_s.second;
        lines.erase(std::remove_if(lines.begin(), lines.end(),
                                   boost::bind(bad_line,
                                               boost::ref(common),
                                               _1)),
                    lines.end());
        Strings lines_2 = remove_empty_sections(lines);
        lines.swap(lines_2);
    }
}

typedef std::map<std::string, std::string> Name2Key;

static std::string link_names(std::string text,
                              const Name2Key& n2k) {
    BOOST_FOREACH (const Name2Key::value_type& n_k, n2k) {
        const std::string& name = n_k.first;
        const std::string& key = n_k.second;
        std::string l = "<a href='#" + key + "'>" + name + "</a>";
        using namespace boost::algorithm;
        replace_all(text, " " + name + ":", " " + l + ":");
        replace_all(text, "\n" + name + ":", "\n" + l + ":");
    }
    return text;
}

void html_all_processors(Meta* m, std::string out_fname) {
    Meta& meta = *m;
    typedef boost::shared_ptr<std::ostream> OPtr;
    OPtr out = name_to_ostream(out_fname);
    std::ostream& o = *out;
    std::string n = "\n";
    o << "<table border='1'>" << n;
    o << "<tr>" << n;
    o << "<td>Summary</td>" << n;
    o << "<td>Help</td>" << n;
    o << "</tr>" << n;
    Key2Strings k2s;
    Name2Key n2k;
    BOOST_FOREACH (std::string key, meta.keys()) {
        SharedProcessor p = meta.get(key);
        add_key(k2s, p.get());
        n2k[p->key()] = p->key();
        n2k[p->name()] = p->key();
    }
    remove_common(k2s);
    BOOST_FOREACH (std::string key, meta.keys()) {
        SharedProcessor p = meta.get(key);
        o << "<tr valign='top'>" << n;
        o << "<td>" << n;
        o << "<a name='" << key << "'></a>" << n;
        o << "<b>" << p->key() << "</b>" << n;
        if (p->key() != p->name()) {
            o << "<br/>" << p->name() << n;
        }
        o << "<pre>";
        set_sstream("tree");
        print_processor_tree("tree", p.get());
        std::string tree = read_file("tree");
        remove_stream("tree");
        using namespace boost::algorithm;
        if (trim_copy(tree).find('\n') != std::string::npos) {
            o << link_names(tree, n2k);
        }
        o << "</pre>" << n;
        print_block_sets(o, p.get());
        o << "</td>" << n;
        o << "<td><pre>" << n;
        const Strings& help = k2s[p->key()];
        o << link_names(boost::join(help, "\n"), n2k) << n;
        o << "</pre></td>" << n;
        o << "</tr>" << n;
    }
    o << "</table>" << n;
}

void html_all_global_options(Meta* m,
                             std::string out_fname) {
    Meta& meta = *m;
    typedef boost::shared_ptr<std::ostream> OPtr;
    OPtr out = name_to_ostream(out_fname);
    std::ostream& o = *out;
    std::string n = "\n";
    o << "<table border='1'>" << n;
    o << "<tr>" << n;
    o << "<td>Global option</td>" << n;
    o << "<td>Value</td>" << n;
    o << "<td>Description</td>" << n;
    o << "</tr>" << n;
    BOOST_FOREACH (std::string opt_name, meta.opts()) {
        std::string opt_value = meta.get_opt(opt_name).to_s();
        std::string opt_d = meta.get_description(opt_name);
        o << "<tr>" << n;
        o << "<td>" << opt_name << "</td>" << n;
        o << "<td>" << opt_value << "</td>" << n;
        o << "<td>" << opt_d << "</td>" << n;
        o << "</tr>" << n;
    }
    o << "</table>" << n;
}

AllProcessors::AllProcessors():
    out_(this, "out", "Output file", true) {
    set_opt_value("out", std::string(":cout"));
}

void AllProcessors::run_impl() const {
    std::string out = opt_value("out").as<std::string>();
    html_all_processors(meta(), out);
}

const char* AllProcessors::name_impl() const {
    return "Print HTML help about all processors";
}

AllOptions::AllOptions():
    out_(this, "out", "Output file", true) {
    set_opt_value("out", std::string(":cout"));
}

void AllOptions::run_impl() const {
    std::string out = opt_value("out").as<std::string>();
    html_all_global_options(meta(), out);
}

const char* AllOptions::name_impl() const {
    return "Print HTML help about all global options";
}

}

