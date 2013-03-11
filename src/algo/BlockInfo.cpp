/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "BlockInfo.hpp"
#include "Block.hpp"

namespace bloomrepeats {

void BlockInfo::print_header(std::ostream& o) const {
    o << "block" << ' ';
    o << "fragments" << ' ';
    o << "cols" << ' ';
    o << "ident-nogap" << ' ';
    o << "ident-gap" << ' ';
    o << "noident-nogap" << ' ';
    o << "noident-gap" << ' ';
    o << "pure-gap" << ' ';
    o << "ident" << std::endl;
}

void BlockInfo::print_block(std::ostream& o, Block* block) const {
    o << block->name() << ' ';
    o << block->size() << ' ';
    AlignmentStat stat;
    block->make_stat(stat);
    o << stat.total << ' ';
    o << stat.ident_nogap << ' ';
    o << stat.ident_gap << ' ';
    o << stat.noident_nogap << ' ';
    o << stat.noident_gap << ' ';
    o << stat.pure_gap << ' ';
    float ident = stat.total ? float(stat.ident_nogap) / float(stat.total) : 0;
    o << ident << std::endl;
}

const char* BlockInfo::name_impl() const {
    return "Print information about blocks";
}

}

