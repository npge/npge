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

UniqueNames::UniqueNames() {
    declare_bs("target", "Target blockset");
}

void UniqueNames::run_impl() const {
    int genomes = genomes_number(*block_set());
    std::set<std::string> names;
    std::string null_name = Block().name(); // 0000 0000
    typedef std::map<std::string, int> String2Int;
    String2Int last_n;
    BOOST_FOREACH (Block* b, *block_set()) {
        if (b->name() == null_name || b->name().empty()) {
            b->set_name(block_name(b, genomes));
        }
        if (names.find(b->name()) != names.end()) {
            std::string orig_name = b->name();
            std::string base_name = orig_name + "_";
            int& i = last_n[orig_name];
            do {
                i += 1;
                b->set_name(base_name + TO_S(i));
            } while (names.find(b->name()) != names.end());
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

const char* UniqueNames::name_impl() const {
    return "Set unique names to all blocks of this block set";
}

}

