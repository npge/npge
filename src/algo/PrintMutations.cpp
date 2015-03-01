/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>

#include "PrintMutations.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

PrintMutations::PrintMutations() {
    declare_bs("target", "Target blockset");
    add_opt("compress", "Compress output by printing . "
            "instead of repeated block or fragment name",
            true);
}

void PrintMutations::find_mutations(const Block* block,
                                    const MutationHandler& func) const {
    TimeIncrementer ti(this);
    std::string cons = block->consensus_string();
    int size = block->size();
    if (size == 1) {
        ASSERT_EQ(block->alignment_length(), cons.size());
    }
    BOOST_FOREACH (Fragment* f, *block) {
        int gaps = 0;
        for (int pos = 0; pos < cons.size(); pos++) {
            char x = f->alignment_at(pos);
            if (size == 1) {
                ASSERT_EQ(x, cons[pos]);
            }
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
    if (m.stop > m.start) {
        ASSERT_EQ(m.change, '-');
        for (int i = m.start; i <= m.stop; i++) {
            o << block->name() << '\t' << f->id() << '\t';
            o << i << '\t' << m.change << '\n';
        }
    } else {
        o << block->name() << '\t' << f->id() << '\t';
        o << m.start << '\t' << m.change << '\n';
    }
}

struct CompressedPrinter {
    std::ostream& o_;
    const Fragment* prev_f_;
    const Block* prev_b_;

    CompressedPrinter(std::ostream& o):
        o_(o), prev_f_(0), prev_b_(0) {
    }

    void print(const Mutation& m) {
        const Fragment* f = m.fragment;
        const Block* block = f->block();
        if (block != prev_b_) {
            o_ << block->name();
            prev_b_ = block;
        } else {
            o_ << '.';
        }
        o_ << '\t';
        if (f != prev_f_) {
            o_ << f->id();
            prev_f_ = f;
        } else {
            o_ << '.';
        }
        o_ << '\t';
        o_ << m.start << '\t';
        if (m.stop > m.start) {
            ASSERT_EQ(m.change, '-');
            o_ << m.stop;
        } else {
            o_ << m.change;
        }
        o_ << '\n';
    }
};

void PrintMutations::print_block(std::ostream& o, Block* block) const {
    TimeIncrementer ti(this);
    if (opt_value("compress").as<bool>()) {
        CompressedPrinter p(o);
        find_mutations(block, boost::bind(
                           &CompressedPrinter::print, &p, _1));
    } else {
        find_mutations(block, boost::bind(print_change,
                                          boost::ref(o), _1));
    }
}

void PrintMutations::print_header(std::ostream& o) const {
    o << "block" << '\t' << "fragment"
      << '\t' << "start" << '\t'
      << "stop(gaps)/change" << '\n';
}

const char* PrintMutations::name_impl() const {
    return "Find all mutations in block";
}

}

