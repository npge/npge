/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "IsPangenome.hpp"
#include "SizeLimits.hpp"
#include "AreBlocksGood.hpp"
#include "UniqueNames.hpp"
#include "Rest.hpp"
#include "AddBlastBlocks.hpp"
#include "Align.hpp"
#include "TrySmth.hpp"
#include "Filter.hpp"
#include "BlockSet.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "Union.hpp"
#include "Subtract.hpp"
#include "block_stat.hpp"
#include "Decimal.hpp"
#include "boundaries.hpp"
#include "hit.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"
#include "global.hpp"

namespace npge {

IsPangenome::IsPangenome() {
    are_blocks_good_ = new AreBlocksGood;
    are_blocks_good_->set_parent(this);
    are_blocks_good_->point_bs("target=target", this);
    align_ = new Align;
    align_->set_parent(this);
    abb_ = new AddBlastBlocks;
    abb_->set_parent(this);
    abb_->point_bs("target=blast-hits", this);
    abb_->point_bs("other=target", this);
    try_join_ = new TrySmth;
    try_join_->set_parent(this);
    try_join_->set_options("--smth-processor:=Joiner");
    try_join_->point_bs("target=joined", this);
    ASSERT_NE(try_join_->block_set(), block_set());
    declare_bs("target", "Target blockset to be tested");
    declare_bs("blast-hits", "Good blast hits");
    declare_bs("all-blast-hits", "All blast hits");
    declare_bs("non-internal-hits", "Non-internal blast hits");
    declare_bs("joined", "Results of joining neighbour blocks");
}

static void remove_non_internal_hits(const BlockSetPtr& hits,
                                     const BlockSetPtr& block_set) {
    SetFc s2f;
    s2f.add_bs(*block_set);
    Blocks hits_blocks(hits->begin(), hits->end());
    BOOST_FOREACH (Block* hit, hits_blocks) {
        if (!is_internal_hit(s2f, hit)) {
            hits->erase(hit);
        }
    }
}

static void fix_self_overlaps_in_hits(const BlockSetPtr& hits) {
    BOOST_FOREACH (Block* hit, *hits) {
        fix_self_overlaps(hit);
    }
}

void IsPangenome::run_impl() const {
    ASSERT_EQ(are_blocks_good_->block_set(), block_set());
    bool good = are_blocks_good_->are_blocks_good();
    //
    UniqueNames un;
    std::ostream& out = are_blocks_good_->get_out();
    try_join_->block_set()->clear();
    Union u;
    u.set_bs("other", block_set());
    u.set_bs("target", try_join_->block_set());
    u.run();
    LiteFilter f;
    f.set_opt_value("min-block", 2);
    f.set_opt_value("min-fragment", 0);
    f.apply(try_join_->block_set());
    try_join_->run();
    Subtract subtract;
    subtract.set_other(block_set());
    subtract.set_block_set(try_join_->block_set());
    subtract.set_opt_value("subtract-equal", true);
    subtract.run();
    f.apply(try_join_->block_set());
    if (!try_join_->block_set()->empty()) {
        good = false;
        out << "Some blocks can be joined" << "\n";
        un.apply(try_join_->block_set());
    }
    //
    abb_->run();
    BlockSetPtr hits = abb_->block_set();
    Union all_hits(hits);
    all_hits.apply(get_bs("all-blast-hits"));
    un.apply(get_bs("all-blast-hits"));
    if (!hits->empty()) {
        align_->apply(hits);
        fix_self_overlaps_in_hits(hits);
        align_->apply(hits);
        Union non_internal_hits(hits);
        non_internal_hits.apply(get_bs("non-internal-hits"));
        un.apply(get_bs("non-internal-hits"));
        if (!hits->empty()) {
            out << "There are " << hits->size() <<
                " non internal hits" << std::endl;
        }
        remove_non_internal_hits(hits, block_set());
        if (!hits->empty()) {
            good = false;
            Boundaries lengths;
            Boundaries sizes;
            Decimals identities;
            BOOST_FOREACH (Block* b, *hits) {
                lengths.push_back(b->alignment_length());
                sizes.push_back(b->size());
                AlignmentStat al_stat;
                make_stat(al_stat, b);
                Decimal identity = block_identity(al_stat);
                identities.push_back(identity);
            }
            double avg_hit_length = avg_element_double(lengths);
            double avg_hit_size = avg_element_double(sizes);
            Decimal avg_hit_identity = avg_element_double(identities);
            out << "There are " << hits->size() << " blast hits "
                << "found on consensuses of blocks.\n"
                << "Average length: " << avg_hit_length << " np.\n"
                << "Average size: " << avg_hit_size << " fragments\n"
                << "Average identity (mapped to orig. blocks): "
                << avg_hit_identity << "\n"
               ;
            un.apply(hits);
        }
    }
    if (good) {
        out << "[good pangenome]" << std::endl;
    } else {
        out << "[not good pangenome]" << std::endl;
    }
}

const char* IsPangenome::name_impl() const {
    return "Print if blockset is good pangenome";
}

}

