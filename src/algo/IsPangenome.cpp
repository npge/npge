/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/join.hpp>

#include "IsPangenome.hpp"
#include "SizeLimits.hpp"
#include "Connector.hpp"
#include "Rest.hpp"
#include "AddBlastBlocks.hpp"
#include "Align.hpp"
#include "MoveGaps.hpp"
#include "CutGaps.hpp"
#include "TrySmth.hpp"
#include "Filter.hpp"
#include "BlockSet.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "Union.hpp"
#include "UniqueNames.hpp"
#include "block_stat.hpp"
#include "boundaries.hpp"
#include "hit.hpp"
#include "block_hash.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

IsPangenome::IsPangenome():
    file_writer_(this, "out-is-pangenome", "Output file with verdict") {
    add_size_limits_options(this);
    move_gaps_ = new MoveGaps();
    move_gaps_->set_parent(this);
    cut_gaps_ = new CutGaps();
    cut_gaps_->set_parent(this);
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
    BOOST_ASSERT(try_join_->block_set() != block_set());
}

static void remove_non_internal_hits(const BlockSetPtr& hits,
                                     const BlockSetPtr& block_set) {
    S2F s2f;
    s2f.add_bs(*block_set);
    typedef std::vector<Block*> Blocks;
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
    bool good = true;
    UniqueNames un;
    Connector c;
    c.apply(block_set());
    Rest r(block_set());
    r.run();
    std::ostream& out = file_writer_.output();
    if (!r.block_set()->empty()) {
        good = false;
        out << "Sequences must be covered entirely by blocks. ";
        out << "There are " << r.block_set()->size()
            << " uncovered regions." << std::endl;
    }
    std::vector<std::string> alignmentless_blocks;
    std::vector<std::string> bad_identity_blocks;
    std::vector<std::string> bad_length_blocks;
    std::vector<std::string> bad_cut_gaps_blocks;
    std::vector<std::string> bad_move_gaps_blocks;
    std::vector<std::string> overlaps_blocks;
    std::vector<std::string> self_overlaps_blocks;
    std::vector<std::string> neighbour_unique;
    int min_fragment_length = opt_value("min-fragment").as<int>();
    double min_identity = opt_value("min-identity").as<double>();
    BOOST_FOREACH (Block* b, *block_set()) {
        AlignmentStat al_stat;
        make_stat(al_stat, b);
        if (al_stat.min_fragment_length() < min_fragment_length &&
                b->size() > 1) {
            bad_length_blocks.push_back(b->name());
        }
        if (al_stat.overlapping_fragments()) {
            overlaps_blocks.push_back(b->name());
            if (has_self_overlaps(b)) {
                self_overlaps_blocks.push_back(b->name());
            }
        }
        if (b->size() != 1) {
            if (al_stat.alignment_rows() != b->size()) {
                alignmentless_blocks.push_back(b->name());
            } else {
                double identity = block_identity(al_stat);
                if (identity < min_identity) {
                    bad_identity_blocks.push_back(b->name());
                }
                boost::shared_ptr<Block> copy(b->clone());
                if (move_gaps_->move_gaps(copy.get())) {
                    bad_move_gaps_blocks.push_back(b->name());
                }
                boost::shared_ptr<Block> copy2(b->clone());
                if (cut_gaps_->cut_gaps(copy2.get())) {
                    bad_cut_gaps_blocks.push_back(b->name());
                }
            }
        } else {
            const Fragment* f = b->front();
            for (int ori = -1; ori <= 1; ori += 2) {
                const Fragment* neighbour = f->neighbor(ori);
                if (neighbour && neighbour->block() &&
                        neighbour->block()->size() == 1) {
                    neighbour_unique.push_back(f->id());
                }
            }
        }
    }
    if (!alignmentless_blocks.empty()) {
        good = false;
        out << "Following blocks do not have alignment: "
            << boost::algorithm::join(alignmentless_blocks, " ")
            << ".\n\n";
    }
    if (!bad_identity_blocks.empty()) {
        good = false;
        out << "Following blocks have identity less then "
            << min_identity << ": "
            << boost::algorithm::join(bad_identity_blocks, " ")
            << ".\n\n";
    }
    if (!bad_length_blocks.empty()) {
        good = false;
        out << "Following blocks have fragments with length less then "
            << min_fragment_length << ": "
            << boost::algorithm::join(bad_length_blocks, " ")
            << ".\n\n";
    }
    if (!bad_move_gaps_blocks.empty()) {
        good = false;
        out << "Following blocks have short 'tails' in alignment: "
            << boost::algorithm::join(bad_move_gaps_blocks, " ")
            << ".\n\n";
    }
    if (!bad_cut_gaps_blocks.empty()) {
        good = false;
        out << "Following blocks have end gaps in alignment: "
            << boost::algorithm::join(bad_cut_gaps_blocks, " ")
            << ".\n\n";
    }
    if (!neighbour_unique.empty()) {
        good = false;
        out << "Following unique fragments have unique neighbours: "
            << boost::algorithm::join(neighbour_unique, " ")
            << ".\n\n";
    }
    if (!overlaps_blocks.empty()) {
        good = false;
        out << "Following blocks have fragments overlapping neighbours: "
            << boost::algorithm::join(overlaps_blocks, " ")
            << ".\n\n";
    }
    if (!self_overlaps_blocks.empty()) {
        good = false;
        out << "Following blocks have self overlapping fragments: "
            << boost::algorithm::join(self_overlaps_blocks, " ")
            << ".\n\n";
    }
    //
    try_join_->block_set()->clear();
    Union u;
    u.set_bs("other", block_set());
    u.set_bs("target", try_join_->block_set());
    u.run();
    Filter f(0, 2);
    f.apply(try_join_->block_set());
    size_t hash_1 = blockset_hash(*try_join_->block_set(), workers());
    try_join_->run();
    Rest r2;
    r2.set_bs("other", try_join_->block_set());
    r2.set_bs("target", try_join_->block_set());
    r2.run();
    un.apply(try_join_->block_set());
    c.apply(try_join_->block_set());
    size_t hash_2 = blockset_hash(*try_join_->block_set(), workers());
    if (hash_1 != hash_2) {
        good = false;
        out << "Some blocks can be joined" << "\n";
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
            Floats identities;
            BOOST_FOREACH (Block* b, *hits) {
                lengths.push_back(b->alignment_length());
                sizes.push_back(b->size());
                AlignmentStat al_stat;
                make_stat(al_stat, b);
                double identity = block_identity(al_stat);
                identities.push_back(identity);
            }
            double avg_hit_length = avg_element_double(lengths);
            double avg_hit_size = avg_element_double(sizes);
            double avg_hit_identity = avg_element_double(identities);
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
    return "Print if block set is good pangenome";
}

}

