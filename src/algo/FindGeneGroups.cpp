/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/thread/tss.hpp>

#include "FindGeneGroups.hpp"
#include "FragmentCollection.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "convert_position.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

struct FindGeneGroups::Impl {
    typedef std::vector<Fragment*> Fragments;
    typedef FragmentCollection<Fragment*, Fragments> FC;
    FC fc_;
    typedef std::vector<Block*> Blocks;
    boost::thread_specific_ptr<Blocks> thread_blocks_;
    Blocks blocks_;
    boost::mutex blocks_mutex_;
};

FindGeneGroups::FindGeneGroups():
    BlocksJobs("pangenome") {
    impl_ = new Impl;
}

FindGeneGroups::~FindGeneGroups() {
    delete impl_;
    impl_ = 0;
}

void FindGeneGroups::change_blocks_impl(std::vector<Block*>& ) const {
    impl_->fc_.clear();
    impl_->fc_.add_bs(*get_bs("genes"));
    impl_->fc_.prepare();
}

bool FindGeneGroups::initialize_thread_impl() const {
    impl_->thread_blocks_.reset(new Impl::Blocks);
    return false;
}

typedef std::pair<int, int> Coords; // in pangenome block, min and max pos.
typedef std::map<Fragment*, Coords> F2C;

struct GenesFragmentComp {
    GenesFragmentComp(F2C* f2c):
        f2c_(f2c)
    { }

    bool operator()(Fragment* a, Fragment* b) const {
        return (*f2c_)[a] < (*f2c_)[b];
    }

    F2C* f2c_;
};

bool FindGeneGroups::apply_to_block_impl(Block* block) const {
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
        BOOST_ASSERT(pangenome_fragment);
        int sequence_begin = gene_part->begin_pos();
        int sequence_last = gene_part->last_pos();
        int fr_begin = seq_to_frag(pangenome_fragment, sequence_begin);
        int fr_last = seq_to_frag(pangenome_fragment, sequence_last);
        int pangenome_begin = block_pos(pangenome_fragment,
                fr_begin, block_length);
        int pangenome_last = block_pos(pangenome_fragment,
                fr_last, block_length);
        f2c[gene_part] = std::make_pair(pangenome_begin, pangenome_last);
    }
    std::sort(gene_parts.begin(), gene_parts.end(), GenesFragmentComp(&f2c));
    Fragment* prev = 0;
    Block* gene_group = 0;
    Impl::Blocks& thread_blocks_ = *impl_->thread_blocks_;
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
            gene_group->set_name(block->name() + "_" +
                                 boost::lexical_cast<std::string>(number));
            thread_blocks_.push_back(gene_group);
            prev = 0;
        }
        gene_group->insert(gene_part);
        prev = gene_part;
    }
    return false;
}

bool FindGeneGroups::finish_thread_impl() const {
    boost::mutex::scoped_lock lock(impl_->blocks_mutex_);
    BOOST_FOREACH (Block* b, *impl_->thread_blocks_) {
        impl_->blocks_.push_back(b);
    }
    impl_->thread_blocks_->clear();
    return false;
}

bool FindGeneGroups::finish_work_impl() const {
    bool result = !impl_->blocks_.empty();
    BlockSet& target = *get_bs("target");
    BOOST_FOREACH (Block* b, impl_->blocks_) {
        target.insert(b);
    }
    impl_->blocks_.clear();
    return result;
}

const char* FindGeneGroups::name_impl() const {
    return "Find groups of gene parts according to pangenome";
}

}

