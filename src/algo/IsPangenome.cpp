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
#include "Connector.hpp"
#include "Rest.hpp"
#include "AddBlastBlocks.hpp"
#include "ExternalAligner.hpp"
#include "Filter.hpp"
#include "BlockSet.hpp"
#include "Block.hpp"
#include "block_stat.hpp"
#include "boundaries.hpp"
#include "process.hpp"

namespace bloomrepeats {

void IsPangenome::add_options_impl(po::options_description& desc) const {
    SizeLimits::add_options_impl(desc);
    add_unique_options(desc)
    ("out-is-pangenome", po::value<std::string>(), "Output file with verdict")
   ;
}

void IsPangenome::apply_options_impl(const po::variables_map& vm) {
    SizeLimits::apply_options_impl(vm);
    if (vm.count("out-is-pangenome")) {
        set_output_file(vm["out-is-pangenome"].as<std::string>());
    }
}

bool IsPangenome::run_impl() const {
    bool good = true;
    Connector c;
    c.apply(block_set());
    Rest r(block_set());
    r.run();
    if (!r.block_set()->empty()) {
        good = false;
        output() << "Sequences must be covered entirely by blocks. ";
        output() << "There are " << r.block_set()->size()
                 << " uncovered regions." << std::endl;
    }
    std::vector<std::string> bad_identity_blocks;
    std::vector<std::string> bad_length_blocks;
    std::vector<std::string> overlaps_blocks;
    BOOST_FOREACH (Block* b, *block_set()) {
        AlignmentStat al_stat;
        make_stat(al_stat, b);
        float identity = block_identity(al_stat);
        if (identity < min_identity()) {
            bad_identity_blocks.push_back(b->name());
        }
        if (al_stat.min_fragment_length() < min_fragment_length() &&
                b->size() > 1) {
            bad_length_blocks.push_back(b->name());
        }
        if (al_stat.overlapping_fragments()) {
            overlaps_blocks.push_back(b->name());
        }
    }
    if (!bad_identity_blocks.empty()) {
        good = false;
        output() << "Following blocks have identity less then "
                 << min_identity() << ": "
                 << boost::algorithm::join(bad_identity_blocks, " ")
                 << ".\n\n";
    }
    if (!bad_length_blocks.empty()) {
        good = false;
        output() << "Following blocks have fragments with length less then "
                 << min_fragment_length() << ": "
                 << boost::algorithm::join(bad_length_blocks, " ")
                 << ".\n\n";
    }
    if (!overlaps_blocks.empty()) {
        good = false;
        output() << "Following blocks have fragments overlapping neighbours: "
                 << boost::algorithm::join(overlaps_blocks, " ")
                 << ".\n\n";
    }
    po::options_description desc;
    add_options(desc);
    AddBlastBlocks abb(block_set());
    copy_processor_options(abb, *this);
    int ll = min_fragment_length() * 2;
    abb.set_options("--blast-min-length=" +
                    boost::lexical_cast<std::string>(ll));
    float li = 1 - (1 - min_identity()) / 2;
    abb.set_options("--blast-min-ident=" +
                    boost::lexical_cast<std::string>(li));
    abb.run();
    if (!abb.block_set()->empty()) {
        BlockSetPtr hits = abb.block_set();
        ExternalAligner ea;
        ea.apply(hits);
        Filter f(min_fragment_length());
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
            output() << "There are " << hits->size() << " blast hits "
                     << "found on consensuses of blocks.\n"
                     << "Average length: " << avg_hit_length << " np.\n"
                     << "Average size: " << avg_hit_size << " fragments\n"
                     << "Average identity (mapped to orig. blocks): "
                     << avg_hit_identity << "\n"
                    ;
        }
    }
    if (good) {
        output() << "[good pangenome]" << std::endl;
    } else {
        output() << "[not good pangenome]" << std::endl;
    }
    return true; // Connector
}

const char* IsPangenome::name_impl() const {
    return "Print if block set is good pangenome";
}

}

