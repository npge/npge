/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <iostream>
#include <vector>
#include <boost/foreach.hpp>

#include "Meta.hpp"
#include "Processor.hpp"
#include "process.hpp"
#include "global.hpp"

using namespace bloomrepeats;

int main() {
    std::ostream& o = std::cout;
    std::string n = "\n";
    o << "<table border='1'>" << n;
    o << "<tr>" << n;
    o << "<td>Key</td>" << n;
    o << "<td>Name</td>" << n;
    o << "<td>Block sets</td>" << n;
    o << "<td>Help</td>" << n;
    o << "<td>Tree</td>" << n;
    o << "</tr>" << n;
    Meta meta;
    BOOST_FOREACH (std::string key, meta.keys()) {
        SharedProcessor p = meta.get(key);
        o << "<tr>" << n;
        o << "<td>" << p->key() << "</td>" << n;
        o << "<td>" << p->name() << "</td>" << n;
        o << "<td>TODO</td>" << n;
        o << "<td><pre>";
        print_help(":cout", p.get());
        o << "</pre></td>" << n;
        o << "<td><pre>";
        print_processor_tree(":cout", p.get());
        o << "</pre></td>" << n;
        o << "</tr>" << n;

    }
}

