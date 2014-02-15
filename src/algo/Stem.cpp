/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <set>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/cast.hpp>

#include "Stem.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

Stem::Stem():
        exact_(false) {
    add_opt("exact", "Require exactly one fragment in each genome", false);
}

bool Stem::is_good_block(const Block* block) const {
    std::set<std::string> genomes;
    BOOST_FOREACH (const Fragment* f, *block) {
        std::string genome = f->seq()->genome();
        if (exact_ && genomes.find(genome) != genomes.end()) {
            return false;
        }
        genomes.insert(genome);
    }
    BOOST_FOREACH (const std::string& genome, genomes_) {
        if (genomes.find(genome) == genomes.end()) {
            return false;
        }
    }
    return true;
}

void Stem::calculate_genomes() const {
    genomes_.clear();
    BOOST_FOREACH (const SequencePtr& seq, block_set()->seqs()) {
        BOOST_ASSERT(!seq->genome().empty());
        genomes_.push_back(seq->genome());
    }
    genomes_.sort_unique();
}

bool Stem::initialize_work_impl() const {
    calculate_genomes();
    exact_ = opt_value("exact").as<bool>();
}

class StemData : public ThreadData {
public:
    std::vector<Block*> blocks_to_erase;
};

ThreadData* Stem::before_thread_impl() const {
    return new StemData;
}

bool Stem::process_block_impl(Block* block, ThreadData* d) const {
    if (!is_good_block(block)) {
        StemData* data = boost::polymorphic_downcast<StemData*>(d);
        data->blocks_to_erase.push_back(block);
    }
    return true; // TODO
}

bool Stem::after_thread_impl(ThreadData* d) const {
    StemData* data = boost::polymorphic_downcast<StemData*>(d);
    BlockSet& target = *block_set();
    BOOST_FOREACH (Block* block, data->blocks_to_erase) {
        target.erase(block);
    }
    return true; // TODO
}

const char* Stem::name_impl() const {
    return "Filter out blocks not represented in at least one of genomes.";
}

}

