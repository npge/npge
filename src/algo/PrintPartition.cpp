/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
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
    o << "sequence\t";
    o << "sequence_start\t";
    o << "sequence_stop\t";
    o << "gene\t";
    o << "gene_block_start\t";
    o << "gene_block_stop\t";
    o << "npg_block\t";
    o << "npg_block_start\t";
    o << "npg_block_stop\t";
    o << std::endl;
}

struct GenesCmp {
    int ori_;

    GenesCmp(int ori):
        ori_(ori) {
    }

    bool operator()(const Fragment* a,
                    const Fragment* b) const {
        if (ori_ == -1) {
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
                      int length_before) {
    Block* npg_block = npg->block();
    int npg_length = npg_block->alignment_length();
    //
    int gene_start = length_before;
    int gene_stop = gene_start + overlap->length() - 1;
    int seq_start = overlap->begin_pos();
    int seq_stop = overlap->last_pos();
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
        // use second word from full gene name
        gene_name = words[1];
    }
    o << overlap->seq()->name() << '\t';
    o << seq_start << '\t';
    o << seq_stop << '\t';
    o << gene_name << '\t';
    o << gene_start << '\t';
    o << gene_stop << '\t';
    o << npg_block->name() << '\t';
    o << npg_block_start << '\t';
    o << npg_block_stop << '\t';
    o << std::endl;
}

void printGenePart(std::ostream& o, Fragment* gene_part,
        int length_before, const VectorFc& fc) {
    Block* gene = gene_part->block();
    std::vector<Fragment> overlaps;
    fc.find_overlaps(overlaps, gene_part);
    std::sort(overlaps.begin(), overlaps.end(),
            GenesCmp(gene_part->ori()));
    BOOST_FOREACH (Fragment& overlap, overlaps) {
        std::vector<Fragment*> npg_vec;
        fc.find_overlap_fragments(npg_vec, &overlap);
        ASSERT_EQ(npg_vec.size(), 1);
        Fragment* npg = npg_vec[0];
        printGeneSubPart(o, &overlap, npg, gene,
                         length_before);
        length_before += overlap.length();
    }
}

void PrintPartition::print_block(std::ostream& o,
                                 Block* gene) const {
    Fragments gene_parts(gene->begin(), gene->end());
    std::sort(gene_parts.begin(), gene_parts.end(),
            GenesCmp(gene->front()->ori()));
    //
    int length_before = 0;
    BOOST_FOREACH (Fragment* gene_part, gene_parts) {
        printGenePart(o, gene_part, length_before,
                      impl_->fc_);
        length_before += gene_part->length();
    }
}

const char* PrintPartition::name_impl() const {
    return "Print overlaps between two blocksets as table";
}

}

