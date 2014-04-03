/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
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

#include "Meta.hpp"
#include "Processor.hpp"
#include "process.hpp"
#include "name_to_stream.hpp"
#include "read_file.hpp"
#include "global.hpp"
#include "throw_assert.hpp"

using namespace bloomrepeats;

void print_block_sets(const std::string& output,
                      Processor* p) {
    boost::shared_ptr<std::ostream> output_ptr;
    output_ptr = name_to_ostream(output);
    std::ostream& out = *output_ptr;
    Strings block_sets;
    p->get_block_sets(block_sets);
    out << "<ul>";
    BOOST_FOREACH (const std::string& bs_name, block_sets) {
        std::string descr = p->bs_description(bs_name);
        if (!descr.empty()) {
            out << "<li>";
            out << bs_name << ": " << descr;
            out << "</li>" << "\n";
        }
    }
    out << "</ul>";
}

typedef std::map<std::string, Strings> Key2Strings;

void add_key(Key2Strings& k2s, Processor* p) {
    std::string key = p->key();
    set_sstream("help");
    print_help("help", p);
    std::string contents = read_file("help");
    remove_stream("help");
    using namespace boost::algorithm;
    Strings& lines = k2s[key];
    split(lines, contents, is_any_of("\n"));
}

std::string replace_multi_spaces(std::string line) {
    int length = -1;
    while (line.length() != length) {
        length = line.length();
        using namespace boost::algorithm;
        replace_all(line, "  ", " ");
    }
    return line;
}

typedef std::set<std::string> StringSet;

bool bad_line(const StringSet& cmn, const std::string& line) {
    return cmn.find(replace_multi_spaces(line)) != cmn.end() ||
           line.find("-timing ") != std::string::npos ||
           line.find("-workers ") != std::string::npos;
}

Strings remove_empty_sections(const Strings& lines) {
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

void remove_common(Key2Strings& k2s) {
    BOOST_ASSERT(!k2s.empty());
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

int main() {
    std::ostream& o = std::cout;
    std::string n = "\n";
    o << "<table border='1'>" << n;
    o << "<tr>" << n;
    o << "<td>Name</td>" << n;
    o << "<td>Block sets</td>" << n;
    o << "<td>Help</td>" << n;
    o << "<td>Tree</td>" << n;
    o << "</tr>" << n;
    Meta meta;
    Key2Strings k2s;
    BOOST_FOREACH (std::string key, meta.keys()) {
        SharedProcessor p = meta.get(key);
        add_key(k2s, p.get());
    }
    remove_common(k2s);
    BOOST_FOREACH (std::string key, meta.keys()) {
        SharedProcessor p = meta.get(key);
        o << "<tr>" << n;
        o << "<td><b>" << p->key() << "</b>" << n;
        if (p->key() != p->name()) {
            o << "<br/>" << p->name() << n;
        }
        o << "</td>" << n;
        o << "<td>";
        print_block_sets(":cout", p.get());
        o << "</td>" << n;
        o << "<td><pre>" << n;
        const Strings& help = k2s[p->key()];
        o << boost::join(help, "\n") << n;
        o << "</pre></td>" << n;
        o << "<td><pre>";
        print_processor_tree(":cout", p.get());
        o << "</pre></td>" << n;
        o << "</tr>" << n;
    }
}

