/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <sstream>
#include <iostream>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "process.hpp"
#include "Meta.hpp"
#include "po.hpp"
#include "name_to_stream.hpp"
#include "global.hpp"

namespace npge {

void print_processor_tree(const std::string& output,
                          Processor* processor, int indent) {
    boost::shared_ptr<std::ostream> output_ptr;
    output_ptr = name_to_ostream(output);
    std::ostream& out = *output_ptr;
    if (!processor->parent()) {
        out << "\n";
    }
    const int SPACES_IN_TAB = 4;
    std::string tab(SPACES_IN_TAB * indent, ' ');
    out << tab << processor->key();
    out << ": " << processor->name() << "\n";
    BOOST_FOREACH (Processor* child, processor->children()) {
        print_processor_tree(output, child, indent + 1);
    }
}

void print_help(const std::string& output, const Processor* processor,
                const std::string& app, const std::string& positional) {
    boost::shared_ptr<std::ostream> output_ptr = name_to_ostream(output);
    std::ostream& out = *output_ptr;
    out << "Usage:" << std::endl;
    out << app << " [options]";
    if (!positional.empty()) {
        using namespace boost::algorithm;
        std::string::const_iterator it = std::find_if(positional.begin(),
                                         positional.end(), !is_any_of("-"));
        if (it != positional.end()) {
            out << ' ' << std::string(it, positional.end());
        }
    }
    po::options_description desc(processor->name());
    add_general_options(desc);
    processor->add_options(desc);
    out << std::endl << std::endl;
    std::stringstream ss;
    ss << desc;
    std::string desc_str = ss.str();
    std::string desc_str1;
    int s = desc_str.size();
    for (int i = 0; i < s - 5; i++) {
        char c = desc_str[i];
        if (c == '-' && desc_str[i + 1] == '-' &&
                desc_str[i + 3] == ' ') {
            desc_str1 += ' ';
        } else {
            desc_str1 += desc_str[i];
        }
    }
    for (int i = s - 5; i < s; i++) {
        desc_str1 += desc_str[i];
    }
    out << desc_str1 << std::endl;
}

static void print_config_of_section(
    std::string section,
    std::ostream& o,
    const Meta* meta) {
    std::string s_name = section.empty() ? "other" : section;
    std::string line = "-- " + s_name + " --";
    std::string decor(line.size(), '-');
    o << decor << "\n" << line << "\n" << decor << "\n\n";
    BOOST_FOREACH (std::string opt_name,
                  meta->opts_of_section(section)) {
        AnyAs value = meta->get_opt(opt_name);
        std::string opt_value = value.to_lua();
        std::string opt_d = meta->get_description(opt_name);
        o << "-- " << opt_d << "\n";
        o << opt_name << " = " << opt_value << "\n\n";
    }
}

void print_config(const std::string& out, const Meta* meta) {
    typedef boost::shared_ptr<std::ostream> OStreamPtr;
    OStreamPtr output_ptr = name_to_ostream(out);
    std::ostream& o = *output_ptr;
    print_config_of_section("main", o, meta);
    BOOST_FOREACH (std::string section, meta->sections()) {
        if (section != "main" && section != "") {
            print_config_of_section(section, o, meta);
        }
    }
    print_config_of_section("", o, meta);
}

}

