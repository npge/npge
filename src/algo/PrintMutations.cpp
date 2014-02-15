/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>

#include "PrintMutations.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "global.hpp"

namespace bloomrepeats {

PrintMutations::PrintMutations()
{ }

void PrintMutations::find_mutations(const Block* block,
                                    const MutationHandler& func) const {
    std::string cons = block->consensus_string();
    int gaps = 0;
    BOOST_FOREACH (Fragment* f, *block) {
        for (int pos = 0; pos < cons.size(); pos++) {
            char x = f->alignment_at(pos);
            if (x == '\0') {
                gaps += 1;
            }
            if (x != '\0' && gaps) {
                Mutation m;
                m.fragment = f;
                m.start = pos - gaps;
                m.stop = pos - 1;
                m.change = '-';
                func(m);
                gaps = 0;
            }
            if (x != '\0' && x != cons[pos]) {
                Mutation m;
                m.fragment = f;
                m.start = pos;
                m.stop = pos;
                m.change = x;
                func(m);
            }
        }
    }
}

static void print_change(std::ostream& o, const Mutation& m) {
    const Fragment* f = m.fragment;
    const Block* block = f->block();
    o << block->name() << '\t' << f->id() << '\t';
    o << m.start << '\t' << m.stop << '\t' << m.change << '\n';
}

void PrintMutations::print_block(std::ostream& o, Block* block) const {
    find_mutations(block, boost::bind(print_change, boost::ref(o), _1));
}

void PrintMutations::print_header(std::ostream& o) const {
    o << "block" << '\t' << "fragment"
      << '\t' << "start" << '\t' << "stop"
      << '\t' << "change" << '\n';
}

}

