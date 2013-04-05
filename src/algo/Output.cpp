/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "Output.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"

namespace bloomrepeats {

Output::Output(const std::string& prefix):
    export_alignment_(true) {
    set_prefix(prefix);
}

void Output::add_options_impl(po::options_description& desc) const {
    AbstractOutput::add_options_impl(desc);
    bloomrepeats::add_unique_options(desc)
    ("export-alignment", po::value<bool>()->default_value(export_alignment()),
     "use alignment information if available")
   ;
}

void Output::apply_options_impl(const po::variables_map& vm) {
    AbstractOutput::apply_options_impl(vm);
    if (vm.count("export-alignment")) {
        set_export_alignment(vm["export-alignment"].as<bool>());
    }
}

static struct FragmentCompareName2 {
    bool operator()(const Fragment* a, const Fragment* b) const {
        typedef boost::tuple<size_t, size_t, int, const std::string&> Tie;
        return Tie(a->min_pos(), a->max_pos(), a->ori(), a->seq()->name()) <
               Tie(b->min_pos(), b->max_pos(), b->ori(), b->seq()->name());
    }
} fcn2;

void Output::print_block(std::ostream& o, Block* block) const {
    std::vector<Fragment*> fragments(block->begin(), block->end());
    std::sort(fragments.begin(), fragments.end(), fcn2);
    BOOST_FOREACH (Fragment* fr, fragments) {
        o << '>';
        fr->print_header(o);
        o << std::endl;
        fr->print_contents(o, export_alignment() ? '-' : 0x00);
        o << std::endl;
    }
    o << std::endl;
}

const char* Output::name_impl() const {
    return "Output block set";
}

}

