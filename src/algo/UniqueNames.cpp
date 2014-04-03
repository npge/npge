/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "UniqueNames.hpp"
#include "Block.hpp"
#include "Sequence.hpp"
#include "BlockSet.hpp"
#include "rand_name.hpp"
#include "block_hash.hpp"
#include "to_s.hpp"

namespace bloomrepeats {

void UniqueNames::run_impl() const {
    typedef std::set<std::string> StringSet;
    StringSet all_genomes;
    BOOST_FOREACH (const SequencePtr& seq, block_set()->seqs()) {
        all_genomes.insert(seq->genome());
    }
    std::set<std::string> names;
    std::string null_name = Block().name(); // 0000 0000
    BOOST_FOREACH (Block* b, *block_set()) {
        if (b->name() == null_name || b->name().empty()) {
            b->set_name(block_name(b, all_genomes.size()));
        }
        if (names.find(b->name()) != names.end()) {
            std::string base_name = b->name() + "_";
            for (int i = 1; names.find(b->name()) != names.end();
                    i++) {
                b->set_name(base_name + TO_S(i));
            }
        }
        names.insert(b->name());
    }
    names.clear();
    BOOST_FOREACH (const SequencePtr& seq, block_set()->seqs()) {
        while (seq->name().empty() ||
                names.find(seq->name()) != names.end()) {
            const int RAND_SEQ_NAME_LENGTH = 8;
            seq->set_name(rand_name(RAND_SEQ_NAME_LENGTH));
        }
        names.insert(seq->name());
    }
}

}

