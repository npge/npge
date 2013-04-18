/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <string>
#include <map>
#include <set>
#include <vector>
#include <boost/foreach.hpp>

#include "Stem.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

typedef std::set<Block*> AffectedBlocks;
typedef std::map<std::string, AffectedBlocks> Genome2AffectedBlocks;

bool Stem::run_impl() const {
    Genome2AffectedBlocks g2ab;
    BOOST_FOREACH (Block* b, *block_set()) {
        BOOST_FOREACH (Fragment* f, *b) {
            BOOST_ASSERT(f->seq());
            std::string genome = f->seq()->genome();
            BOOST_ASSERT(!genome.empty());
            g2ab[genome].insert(b);
        }
    }
    std::set<std::string> genomes;
    BOOST_FOREACH (const SequencePtr& seq, block_set()->seqs()) {
        BOOST_ASSERT(!seq->genome().empty());
        genomes.insert(seq->genome());
    }
    BOOST_FOREACH (const Genome2AffectedBlocks::value_type& seq_and_ab, g2ab) {
        const std::string& genome = seq_and_ab.first;
        BOOST_ASSERT(genomes.find(genome) != genomes.end());
    }
    bool result = false;
    BlockSet& bs = *block_set();
    BOOST_FOREACH (const SequencePtr& seq, block_set()->seqs()) {
        const AffectedBlocks& af = g2ab[seq->genome()];
        const std::vector<Block*> blocks((bs.begin()), bs.end());
        BOOST_FOREACH (Block* b, blocks) {
            if (af.find(b) == af.end()) {
                bs.erase(b);
                result = true;
            }
        }
    }
    return result;
}

const char* Stem::name_impl() const {
    return "Filter out blocks not represented in at least one of genomes.";
}

}

