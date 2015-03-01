/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <algorithm>
#include <boost/foreach.hpp>

#include "FindGeneGroups.hpp"
#include "FragmentCollection.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"
#include "global.hpp"

namespace npge {

struct FindGeneGroups::Impl {
    VectorFc fc_;
    Blocks blocks_;
};

FindGeneGroups::FindGeneGroups():
    BlocksJobs("pangenome") {
    impl_ = new Impl;
    declare_bs("target", "Where gene groups are written to");
    declare_bs("genes", "Blocks of this blockset represent genes");
    declare_bs("pangenome", "Similar parts of genomes");
}

FindGeneGroups::~FindGeneGroups() {
    delete impl_;
    impl_ = 0;
}

void FindGeneGroups::initialize_work_impl() const {
    impl_->fc_.clear();
    impl_->fc_.add_bs(*get_bs("genes"));
    impl_->fc_.prepare();
}

struct GeneGroupsData : public ThreadData {
    Blocks thread_blocks_;
};

ThreadData* FindGeneGroups::before_thread_impl() const {
    return new GeneGroupsData;
}

typedef std::pair<int, int> Coords; // in pangenome block, min and max pos.
typedef std::map<Fragment*, Coords> F2C;

struct GenesFragmentComp {
    GenesFragmentComp(F2C* f2c):
        f2c_(f2c) {
    }

    bool operator()(Fragment* a, Fragment* b) const {
        return (*f2c_)[a] < (*f2c_)[b];
    }

    F2C* f2c_;
};

void FindGeneGroups::process_block_impl(Block* block,
                                        ThreadData* td) const {
    // block from pangenome
    int block_length = block->alignment_length();
    std::vector<Fragment*> gene_parts;
    BOOST_FOREACH (Fragment* fragment, *block) {
        impl_->fc_.find_overlap_fragments(gene_parts, fragment);
    }
    F2C f2c;
    BOOST_FOREACH (Fragment* gene_part, gene_parts) {
        Fragment* pangenome_fragment = 0;
        BOOST_FOREACH (Fragment* f, *block) {
            if (gene_part->is_subfragment_of(*f)) {
                pangenome_fragment = f;
                break;
            }
        }
        ASSERT_TRUE(pangenome_fragment);
        int sequence_begin = gene_part->begin_pos();
        int sequence_last = gene_part->last_pos();
        int fr_begin = seq_to_frag(pangenome_fragment, sequence_begin);
        int fr_last = seq_to_frag(pangenome_fragment, sequence_last);
        int pangenome_begin = block_pos(pangenome_fragment,
                                        fr_begin, block_length);
        int pangenome_last = block_pos(pangenome_fragment,
                                       fr_last, block_length);
        int pangenome_min = std::min(pangenome_begin, pangenome_last);
        int pangenome_max = std::max(pangenome_begin, pangenome_last);
        f2c[gene_part] = std::make_pair(pangenome_min, pangenome_max);
    }
    std::sort(gene_parts.begin(), gene_parts.end(), GenesFragmentComp(&f2c));
    Fragment* prev = 0;
    Block* gene_group = 0;
    GeneGroupsData* ggd;
    ggd = boost::polymorphic_downcast<GeneGroupsData*>(td);
    Blocks& thread_blocks_ = ggd->thread_blocks_;
    int number = 0;
    BOOST_FOREACH (Fragment* gene_part, gene_parts) {
        if (prev && f2c[gene_part].first > f2c[prev].second) {
            // no overlap on block
            gene_group = 0;
        }
        if (!gene_group) {
            number += 1;
            gene_group = new Block;
            gene_group->set_weak(true);
            gene_group->set_name(block->name() + "g" +
                                 TO_S(number));
            thread_blocks_.push_back(gene_group);
            prev = 0;
        }
        gene_group->insert(gene_part);
        prev = gene_part;
    }
}

void FindGeneGroups::after_thread_impl(ThreadData* td) const {
    GeneGroupsData* ggd;
    ggd = boost::polymorphic_downcast<GeneGroupsData*>(td);
    Blocks& thread_blocks_ = ggd->thread_blocks_;
    BlockSet& target = *get_bs("target");
    BOOST_FOREACH (Block* b, thread_blocks_) {
        target.insert(b);
    }
}

const char* FindGeneGroups::name_impl() const {
    return "Find groups of gene parts according to pangenome";
}

}

