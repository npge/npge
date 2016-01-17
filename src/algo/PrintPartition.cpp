/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2016 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <ostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>

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
    VectorFc fc_;
};

PrintPartition::PrintPartition() {
    impl_ = new Impl;
    declare_bs("genes", "Genes");
    declare_bs("npg", "Pangenome");
    set_block_set_name("genes");
    add_opt("group-by-gene",
            "Print all records about one gene on one line",
            false);
}

PrintPartition::~PrintPartition() {
    delete impl_;
    impl_ = 0;
}

void PrintPartition::prepare() const {
    impl_->fc_.clear();
    impl_->fc_.add_bs(*get_bs("npg"));
    impl_->fc_.prepare();
}

void PrintPartition::print_header(std::ostream& o) const {
    bool group_by_gene =
        opt_value("group-by-gene").as<bool>();
    o << "sequence\t";
    o << "sequence_start\t";
    o << "sequence_stop\t";
    o << "gene\t";
    o << "npg_block\t";
    o << "npg_block_min\t";
    o << "npg_block_max\t";
    o << "npg_block_ori\t";
    o << "gene_block_start\t";
    o << "gene_block_stop\t";
    o << std::endl;
}

struct GenesCmp {
    int ori_;

    GenesCmp(int ori):
        ori_(ori) {
    }

    bool operator()(const Fragment* a,
                    const Fragment* b) const {
        if (ori_ == 1) {
            return a->min_pos() < b->min_pos();
        } else {
            return b->min_pos() < a->min_pos();
        }
    }

    bool operator()(const Fragment& a,
                    const Fragment& b) const {
        const GenesCmp& self = *this;
        return self(&a, &b);
    }
};

void printGeneSubPart(std::ostream& o, Fragment* overlap,
                      Fragment* npg, Block* gene,
                      Fragment* gene_part,
                      int length_before,
                      int part_index, int nparts,
                      bool group_by_gene,
                      bool multifragment_gene) {
    ASSERT_EQ(overlap->ori(), gene_part->ori());
    Block* npg_block = npg->block();
    int npg_length = npg_block->alignment_length();
    //
    int gene_start = length_before;
    int gene_stop = gene_start + overlap->length() - 1;
    int seq_start = overlap->begin_pos();
    int seq_stop = overlap->last_pos();
    int printable_seq_start = seq_start;
    int printable_seq_stop = seq_stop;
    if (group_by_gene) {
        printable_seq_start = gene_part->begin_pos();
        printable_seq_stop = gene_part->last_pos();
    }
    int npg_fr_start = seq_to_frag(npg, seq_start);
    int npg_fr_stop = seq_to_frag(npg, seq_stop);
    int npg_block_start = block_pos(npg, npg_fr_start,
                                    npg_length);
    int npg_block_stop = block_pos(npg, npg_fr_stop,
                                   npg_length);
    // print only locus_tag instead of full gene name
    // see https://github.com/npge/npge/issues/23
    std::string gene_name = gene->name();
    Strings words;
    using namespace boost::algorithm;
    split(words, gene_name, isspace, token_compress_on);
    if (words.size() >= 3) {
        // use third word from full gene name
        gene_name = words[2];
    }
    if (!group_by_gene || part_index == 0) {
        o << overlap->seq()->name() << '\t';
        o << printable_seq_start << '\t';
        o << printable_seq_stop << '\t';
        if (multifragment_gene) {
            o << '*';
        }
        o << gene_name << '\t';
    }
    o << npg_block->name() << '\t';
    o << std::min(npg_block_start, npg_block_stop) << '\t';
    o << std::max(npg_block_start, npg_block_stop) << '\t';
    int ori = (npg_block_stop - npg_block_start >= 0) ? 1 : (-1);
    o << ori << '\t';
    o << gene_start << '\t';
    o << gene_stop;
    if (!group_by_gene || part_index == nparts - 1) {
        o << std::endl;
    } else {
        o << '\t';
    }
}

void printGenePart(std::ostream& o, Fragment* gene_part,
        int length_before, const VectorFc& fc,
        bool group_by_gene, bool multifragment_gene) {
    Block* gene = gene_part->block();
    std::vector<Fragment> overlaps;
    fc.find_overlaps(overlaps, gene_part);
    std::sort(overlaps.begin(), overlaps.end(),
            GenesCmp(gene_part->ori()));
    int part_index = 0;
    int nparts = overlaps.size();
    BOOST_FOREACH (Fragment& overlap, overlaps) {
        std::vector<Fragment*> npg_vec;
        fc.find_overlap_fragments(npg_vec, &overlap);
        ASSERT_EQ(npg_vec.size(), 1);
        Fragment* npg = npg_vec[0];
        if (overlap.ori() != gene_part->ori()) {
            overlap.inverse();
        }
        printGeneSubPart(o, &overlap, npg, gene,
                         gene_part,
                         length_before,
                         part_index, nparts,
                         group_by_gene,
                         multifragment_gene);
        length_before += overlap.length();
        part_index += 1;
    }
}

void PrintPartition::print_block(std::ostream& o,
                                 Block* gene) const {
    Fragments gene_parts(gene->begin(), gene->end());
    std::sort(gene_parts.begin(), gene_parts.end(),
            GenesCmp(gene->front()->ori()));
    //
    bool group_by_gene =
        opt_value("group-by-gene").as<bool>();
    int length_before = 0;
    bool multifragment_gene = gene_parts.size() > 1;
    BOOST_FOREACH (Fragment* gene_part, gene_parts) {
        printGenePart(o, gene_part, length_before,
                      impl_->fc_, group_by_gene,
                      multifragment_gene);
        length_before += gene_part->length();
    }
}

const char* PrintPartition::name_impl() const {
    return "Print overlaps between two blocksets as table";
}

}

