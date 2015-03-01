/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>

#include "AllOptions.hpp"
#include "Meta.hpp"
#include "name_to_stream.hpp"
#include "read_file.hpp"
#include "global.hpp"
#include "throw_assert.hpp"

namespace npge {

static void gopts_of_section(std::string section, Meta* m,
                             std::ostream& o) {
    Meta& meta = *m;
    std::string n = "\n";
    if (section == "config") {
        return;
    }
    if (meta.opts_of_section(section).empty()) {
        return;
    }
    o << "<tr>" << n;
    o << "<td colspan='3' align='center' bgcolor='lightgray'>";
    o << (section.empty() ? "other" : section);
    o << "</td></tr>" << n;
    BOOST_FOREACH (std::string opt_name,
                  meta.opts_of_section(section)) {
        AnyAs v = meta.get_opt(opt_name);
        std::string opt_value = v.to_s();
        if (v.type() == typeid(bool)) {
            opt_value = v.to_lua();
        }
        std::string opt_d = meta.get_description(opt_name);
        o << "<tr>" << n;
        o << "<td>" << opt_name << "</td>" << n;
        o << "<td>" << opt_value << "</td>" << n;
        o << "<td>" << opt_d << "</td>" << n;
        o << "</tr>" << n;
    }
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
    gopts_of_section("main", m, o);
    BOOST_FOREACH (std::string section, meta.sections()) {
        if (section != "main" && section != "") {
            gopts_of_section(section, m, o);
        }
    }
    gopts_of_section("", m, o);
    o << "</table>" << n;
}

AllOptions::AllOptions():
    out_(this, "out", "Output file", true) {
    set_opt_value("out", std::string(":stdout"));
}

void AllOptions::run_impl() const {
    std::string out = opt_value("out").as<std::string>();
    html_all_global_options(meta(), out);
}

const char* AllOptions::name_impl() const {
    return "Print HTML help about all global options";
}

}

