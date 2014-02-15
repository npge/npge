/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <boost/foreach.hpp>

#include "PrintMutations.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "global.hpp"

namespace bloomrepeats {

PrintMutations::PrintMutations()
{ }

static void print_change(std::ostream& o, const Fragment* f,
                         int start, int stop, char change) {
    const Block* block = f->block();
    o << block->name() << '\t' << f->id() << '\t';
    o << start << '\t' << stop << '\t' << change << '\n';
}

void PrintMutations::print_block(std::ostream& o, Block* block) const {
    std::string cons = block->consensus_string();
    int gaps = 0;
    BOOST_FOREACH (const Fragment* f, *block) {
        for (int pos = 0; pos < cons.size(); pos++) {
            char x = f->alignment_at(pos);
            if (x == '\0') {
                gaps += 1;
            }
            if (x != '\0' && gaps) {
                print_change(o, f, pos - gaps, pos - 1, '-');
                gaps = 0;
            }
            if (x != '\0' && x != cons[pos]) {
                print_change(o, f, pos, pos, x);
            }
        }
    }
}

void PrintMutations::print_header(std::ostream& o) const {
    o << "block" << '\t' << "fragment"
      << '\t' << "start" << '\t' << "stop"
      << '\t' << "change" << '\n';
}

}

