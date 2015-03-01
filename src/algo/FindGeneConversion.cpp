/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <boost/cast.hpp>
#include <boost/foreach.hpp>

#include "FindGeneConversion.hpp"
#include "FragmentDistance.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "Exception.hpp"
#include "cast.hpp"
#include "global.hpp"

namespace npge {

FindGeneConversion::FindGeneConversion() {
    distance_ = new FragmentDistance;
    distance_->set_parent(this);
    set_block_set_name("other");
    declare_bs("other", "Where gene conversion is searched");
    declare_bs("target", "Where results are added");
}

class BlocksData : public ThreadData {
public:
    Blocks blocks;

    ~BlocksData() {
        BOOST_FOREACH (Block* block, blocks) {
            delete block;
        }
        blocks.clear();
    }
};

ThreadData* FindGeneConversion::before_thread_impl() const {
    return new BlocksData;
}

typedef std::pair<Fragment*, Fragment*> FF;
typedef std::map<FF, double> Distances;
typedef std::set<Fragment*> FragmentsSet;
typedef std::map<std::string, Fragments> Genome2Fragments;

static bool try_add_to_group(Fragment* fr,
                             FragmentsSet& conversion_set,
                             const FragmentsSet& external_set,
                             double& max_internal,
                             double& min_external,
                             Distances& dst) { // not const for []
    double max_internal_local = max_internal;
    double min_external_local = min_external;
    BOOST_FOREACH (Fragment* f, conversion_set) {
        max_internal_local = std::max(max_internal_local, dst[FF(fr, f)]);
    }
    BOOST_FOREACH (Fragment* f, external_set) {
        min_external_local = std::max(min_external_local, dst[FF(fr, f)]);
    }
    if (max_internal_local < min_external_local) {
        max_internal = max_internal_local;
        min_external = min_external_local;
        conversion_set.insert(fr);
        return true;
    } else {
        return false;
    }
}

static void find_gene_conversion(
    Blocks& blocks,
    const Fragments& genome,
    const Fragments& all,
    Distances& distances,
    int& conversion_number) {
    if (genome.size() < 2) {
        return;
    }
    FragmentsSet external_set(all.begin(), all.end());
    BOOST_FOREACH (Fragment* f, genome) {
        external_set.erase(f);
    }
    FragmentsSet genome_set(genome.begin(), genome.end());
    while (!genome_set.empty()) {
        Fragment* reference = *genome_set.begin();
        genome_set.erase(reference);
        FragmentsSet conversion_set;
        conversion_set.insert(reference);
        double max_internal = 0.0, min_external = 1.0;
        Fragments genome_vector(genome_set.begin(), genome_set.end());
        BOOST_FOREACH (Fragment* f, genome_vector) {
            if (try_add_to_group(f, conversion_set, external_set,
                                 max_internal, min_external, distances)) {
                genome_set.erase(f);
            }
        }
        if (conversion_set.size() >= 2) {
            const Block* parent = reference->block();
            Block* new_block = new Block;
            new_block->set_weak(true);
            new_block->set_name(parent->name() + "_" + TO_S(conversion_number));
            conversion_number += 1;
            BOOST_FOREACH (Fragment* f, conversion_set) {
                new_block->insert(f);
            }
            blocks.push_back(new_block);
        }
    }
}

void FindGeneConversion::change_blocks_impl(std::vector<Block*>&) const {
    block_set()->add_sequences(other()->seqs());
}

void FindGeneConversion::process_block_impl(Block* block,
        ThreadData* d) const {
    if (block->size() < 3) {
        return;
    }
    Distances distances;
    Fragments fragments(block->begin(), block->end());
    for (int i = 0; i < fragments.size(); i++) {
        Fragment* f1 = fragments[i];
        for (int j = i + 1; j < fragments.size(); j++) {
            Fragment* f2 = fragments[j];
            double ratio = distance_->fragment_distance(f1, f2).ratio();
            distances[std::make_pair(f1, f2)] = ratio;
            distances[std::make_pair(f2, f1)] = ratio;
        }
    }
    BOOST_FOREACH (Fragment* f, fragments) {
        distances[std::make_pair(f, f)] = 0.0;
    }
    Genome2Fragments genome2fragments;
    BOOST_FOREACH (Fragment* f, fragments) {
        const Sequence* seq = f->seq();
        if (!seq) {
            throw Exception("Fragment without sequence");
        }
        std::string genome = seq->genome();
        if (genome.empty()) {
            throw Exception("Sequence not from genome: " + seq->name());
        }
        genome2fragments[genome].push_back(f);
    }
    BlocksData* data = boost::polymorphic_downcast<BlocksData*>(d);
    Blocks& blocks = data->blocks;
    int conversion_number = 1;
    BOOST_FOREACH (const Genome2Fragments::value_type& g2f,
                  genome2fragments) {
        const std::string& genome = g2f.first;
        const Fragments& genome_fragments = g2f.second;
        find_gene_conversion(blocks, genome_fragments, fragments, distances,
                             conversion_number);
    }
}

void FindGeneConversion::after_thread_impl(ThreadData* d) const {
    BlocksData* data = boost::polymorphic_downcast<BlocksData*>(d);
    BlockSet& target = *block_set();
    BOOST_FOREACH (Block* block, data->blocks) {
        target.insert(block);
    }
    data->blocks.clear();
}

}

