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
#include "Filter.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "Union.hpp"
#include "block_stat.hpp"
#include "boundaries.hpp"
#include "process.hpp"
#include "hit.hpp"

namespace bloomrepeats {

IsPangenome::IsPangenome():
    file_writer_(this, "out-is-pangenome", "Output file with verdict") {
    add_size_limits_options(this);
    move_gaps_ = new MoveGaps();
    move_gaps_->set_parent(this);
    cut_gaps_ = new CutGaps();
    cut_gaps_->set_parent(this);
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

bool IsPangenome::run_impl() const {
    bool good = true;
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
                float identity = block_identity(al_stat);
                if (identity < min_identity) {
                    bad_identity_blocks.push_back(b->name());
                }
                boost::shared_ptr<Block> copy(Union::clone_block(b));
                if (move_gaps_->move_gaps(copy.get())) {
                    bad_move_gaps_blocks.push_back(b->name());
                }
                if (cut_gaps_->cut_gaps(copy.get())) {
                    bad_cut_gaps_blocks.push_back(b->name());
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
    AddBlastBlocks abb(block_set());
    abb.set_bs("target", get_bs("blast-hits"));
    copy_processor_options(abb, *this);
    int ll = min_fragment_length;
    abb.set_options("--blast-min-length=" +
                    boost::lexical_cast<std::string>(ll));
    double li = min_identity;
    abb.set_options("--blast-min-ident=" +
                    boost::lexical_cast<std::string>(li));
    abb.run();
    if (!abb.block_set()->empty()) {
        BlockSetPtr hits = abb.block_set();
        Align ea;
        ea.apply(hits);
        Filter f(min_fragment_length);
        f.apply(hits);
        remove_non_internal_hits(hits, block_set());
        fix_self_overlaps_in_hits(hits);
        f.apply(hits);
        ea.apply(hits);
        f.apply(hits);
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
                float identity = block_identity(al_stat);
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
        }
    }
    if (good) {
        out << "[good pangenome]" << std::endl;
    } else {
        out << "[not good pangenome]" << std::endl;
    }
    return false;
}

const char* IsPangenome::name_impl() const {
    return "Print if block set is good pangenome";
}

}

