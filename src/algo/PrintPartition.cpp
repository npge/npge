/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "PrintPartition.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "FragmentCollection.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

struct PrintPartition::Impl {
    typedef FragmentCollection<Fragment*, Fragments> FC;
    FC fc_;
};

PrintPartition::PrintPartition() {
    impl_ = new Impl;
    declare_bs("target", "Blocks in which overlaps are searched");
    declare_bs("other", "Blocks with which overlaps are searched");
}

PrintPartition::~PrintPartition() {
    delete impl_;
    impl_ = 0;
}

void PrintPartition::prepare() const {
    impl_->fc_.clear();
    impl_->fc_.add_bs(*other());
    impl_->fc_.prepare();
}

void PrintPartition::print_header(std::ostream& o) const {
    o << "sequence\t";
    o << "sequence_begin\t";
    o << "sequence_last\t";
    o << "orig_block\t";
    o << "target_block\t";
    o << "target_block_begin\t";
    o << "target_block_last\t";
    o << "other_block\t";
    o << "other_block_begin\t";
    o << "other_block_last\t";
    o << std::endl;
}

void PrintPartition::print_block(std::ostream& o, Block* target_block) const {
    int target_block_length = target_block->alignment_length();
    BOOST_FOREACH (Fragment* target, *target_block) {
        Block* orig_block = target->block();
        std::vector<Fragment*> overlap_fragments;
        impl_->fc_.find_overlap_fragments(overlap_fragments, target);
        BOOST_FOREACH (Fragment* other, overlap_fragments) {
            Block* other_block = other->block();
            int other_block_length = other_block->alignment_length();
            Fragment overlap = target->common_fragment(*other);
            int sequence_begin = overlap.begin_pos();
            int sequence_last = overlap.last_pos();
            int target_fr_begin = seq_to_frag(target, sequence_begin);
            int target_fr_last = seq_to_frag(target, sequence_last);
            int target_block_begin = block_pos(target, target_fr_begin,
                                               target_block_length);
            int target_block_last = block_pos(target, target_fr_last,
                                              target_block_length);
            int other_fr_begin = seq_to_frag(other, sequence_begin);
            int other_fr_last = seq_to_frag(other, sequence_last);
            int other_block_begin = block_pos(other, other_fr_begin,
                                              other_block_length);
            int other_block_last = block_pos(other, other_fr_last,
                                             other_block_length);
            o << overlap.seq()->name() << '\t';
            o << sequence_begin << '\t';
            o << sequence_last << '\t';
            o << orig_block->name() << '\t';
            o << target_block->name() << '\t';
            o << target_block_begin << '\t';
            o << target_block_last << '\t';
            o << other_block->name() << '\t';
            o << other_block_begin << '\t';
            o << other_block_last << '\t';
            o << std::endl;
        }
    }
}

const char* PrintPartition::name_impl() const {
    return "Print overlaps between two block sets as table";
}

}

