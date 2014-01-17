/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include "BlockInfo.hpp"
#include "Block.hpp"
#include "block_stat.hpp"

namespace bloomrepeats {

BlockInfo::BlockInfo(const std::string& prefix) {
    set_prefix(prefix);
}

void BlockInfo::print_header(std::ostream& o) const {
    o << "block" << '\t';
    o << "fragments" << '\t';
    o << "cols" << '\t';
    o << "ident-nogap" << '\t';
    o << "ident-gap" << '\t';
    o << "noident-nogap" << '\t';
    o << "noident-gap" << '\t';
    o << "pure-gap" << '\t';
    o << "ident" << '\t';
    o << "GC" << std::endl;
}

void BlockInfo::print_block(std::ostream& o, Block* block) const {
    o << block->name() << '\t';
    o << block->size() << '\t';
    AlignmentStat stat;
    make_stat(stat, block);
    o << stat.total() << '\t';
    o << stat.ident_nogap() << '\t';
    o << stat.ident_gap() << '\t';
    o << stat.noident_nogap() << '\t';
    o << stat.noident_gap() << '\t';
    o << stat.pure_gap() << '\t';
    o << block_identity(stat) << '\t';
    o << stat.gc() << std::endl;
}

const char* BlockInfo::name_impl() const {
    return "Print information about blocks";
}

}

