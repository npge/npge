/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <ostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "PrintGeneGroups.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "FragmentCollection.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace bloomrepeats {

struct PrintGeneGroups::Impl {
    typedef FragmentCollection<Fragment*, Fragments> FC;
    FC fc_;
};

PrintGeneGroups::PrintGeneGroups() {
    impl_ = new Impl;
    declare_bs("target", "Gene groups, weak blocks");
    declare_bs("pangenome", "Similar parts of genomes");
}

PrintGeneGroups::~PrintGeneGroups() {
    delete impl_;
    impl_ = 0;
}

void PrintGeneGroups::prepare() const {
    impl_->fc_.clear();
    impl_->fc_.add_bs(*get_bs("pangenome"));
    impl_->fc_.prepare();
}

void PrintGeneGroups::print_header(std::ostream& o) const {
    o << "group_name\t";
    o << "is_good\t";
    o << "group_fragments\t";
    o << "block\t";
    o << "block_fragments\t";
    o << "block_length\t";
    o << "block_ori\t";
    o << "block_first_min\t";
    o << "block_first_max\t";
    o << "block_last_min\t";
    o << "block_last_max\t";
    o << "has_gene_start\t";
    o << "has_gene_stop\t";
    o << "gene_first_min\t";
    o << "gene_first_max\t";
    o << "gene_last_min\t";
    o << "gene_last_max\t";
    o << std::endl;
}

static void pos_in_gene(const Fragment* gene_part,
                        int& gene_min, int& gene_max) {
    bool search_less = gene_part->ori() == 1;
    Block* gene = gene_part->block();
    int before = 0;
    BOOST_FOREACH (Fragment* f, *gene) {
        BOOST_ASSERT(gene_part->length() == 1 || f->length() == 1
                     /* TODO can't detect ori */ ||
                     f->ori() == gene_part->ori());
        if (f != gene_part) {
            bool less = *f < *gene_part;
            if (search_less == less) {
                before += f->length();
            }
        }
    }
    gene_min = before;
    gene_max = gene_min + gene_part->length() - 1;
    BOOST_ASSERT(gene_min <= gene_max);
}

// start_stop: 1 (start), 2 (stop)
static bool is_gene_boundary(const Fragment* gene_part, int start_stop) {
    if (gene_part->length() == 1) {
        // TODO can't detect ori
        return false;
    }
    Block* gene = gene_part->block();
    bool search_less = start_stop == gene_part->ori();
    BOOST_FOREACH (Fragment* f, *gene) {
        if (f->length() == 1) {
            // TODO can't detect ori
            return false;
        }
        BOOST_ASSERT(f->ori() == gene_part->ori());
        bool less = *f < *gene_part;
        if (search_less == less) {
            return false;
        }
    }
    return true;
}

void PrintGeneGroups::print_block(std::ostream& o, Block* group) const {
    if (group->empty()) {
        return;
    }
    std::vector<Fragment*> fragments1;
    impl_->fc_.find_overlap_fragments(fragments1, group->front());
    BOOST_ASSERT(fragments1.size() == 1);
    Block* block = fragments1.front()->block();
    int block_length = block->alignment_length();
    int block_first_min = block_length;
    int block_first_max = 0;
    int block_last_min = block_length;
    int block_last_max = 0;
    int gene_first_min = -1;
    int gene_first_max = -1;
    int gene_last_min = -1;
    int gene_last_max = -1;
    int ori = -2; // 0 means any, -2 in the beginning
    bool has_gene_start = false;
    bool has_gene_stop = false;
    BOOST_FOREACH (Fragment* gene_part, *group) {
        std::vector<Fragment*> fragments2;
        impl_->fc_.find_overlap_fragments(fragments2, gene_part);
        BOOST_ASSERT(fragments2.size() == 1);
        Fragment* pangenome_fragment = fragments2.front();
        BOOST_ASSERT(pangenome_fragment->block() == block);
        BOOST_ASSERT(gene_part->is_subfragment_of(*pangenome_fragment));
        int sequence_begin = gene_part->begin_pos();
        int sequence_last = gene_part->last_pos();
        int fr_begin = seq_to_frag(pangenome_fragment, sequence_begin);
        int fr_last = seq_to_frag(pangenome_fragment, sequence_last);
        int block_begin = block_pos(pangenome_fragment, fr_begin, block_length);
        int block_last = block_pos(pangenome_fragment, fr_last, block_length);
        int block_min = std::min(block_begin, block_last);
        int block_max = std::max(block_begin, block_last);
        int block_ori = (block_min == block_begin) ? 1 : -1;
        block_first_min = std::min(block_first_min, block_min);
        block_first_max = std::max(block_first_max, block_min);
        block_last_min = std::min(block_last_min, block_max);
        block_last_max = std::max(block_last_max, block_max);
        if (ori == -2) {
            ori = block_ori;
        } else if (ori != block_ori) {
            ori = 0;
        }
        has_gene_start |= is_gene_boundary(gene_part, /* start */ 1);
        has_gene_stop |= is_gene_boundary(gene_part, /* stop */ -1);
        int gene_min, gene_max;
        pos_in_gene(gene_part, gene_min, gene_max);
        if (gene_first_min == -1) {
            gene_first_min = gene_min;
            gene_first_max = gene_min;
            gene_last_min = gene_max;
            gene_last_max = gene_max;
        } else {
            gene_first_min = std::min(gene_first_min, gene_min);
            gene_first_max = std::max(gene_first_max, gene_min);
            gene_last_min = std::min(gene_last_min, gene_max);
            gene_last_max = std::max(gene_last_max, gene_max);
        }
    }
    bool is_good =
        group->size() == block->size() &&
        ori != 0 &&
        block_first_min == block_first_max &&
        block_last_min == block_last_max;
    o << group->name() << '\t';
    o << (is_good ? "good" : "bad") << '\t';
    o << group->size() << '\t';
    o << block->name() << '\t';
    o << block->size() << '\t';
    o << block_length << '\t';
    o << ori << '\t';
    o << block_first_min << '\t';
    o << block_first_max << '\t';
    o << block_last_min << '\t';
    o << block_last_max << '\t';
    o << has_gene_start << '\t';
    o << has_gene_stop << '\t';
    o << gene_first_min << '\t';
    o << gene_first_max << '\t';
    o << gene_last_min << '\t';
    o << gene_last_max << '\t';
    o << std::endl;
}

const char* PrintGeneGroups::name_impl() const {
    return "Print information about gene groups";
}

}

