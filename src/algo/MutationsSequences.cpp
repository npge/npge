/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/cast.hpp>
#include <boost/bind.hpp>

#include "MutationsSequences.hpp"
#include "PrintMutations.hpp"
#include "SeqStorage.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

MutationsSequences::MutationsSequences() {
    set_block_set_name("other");
    print_mutations_ = new PrintMutations;
    print_mutations_->set_parent(this);
    add_opt("mutation-distance",
            "Distance to sequence around mutation to keep", 0);
    add_seq_storage_options(this);
}

typedef std::map<std::string, std::string> Genome2Str;
typedef std::map<Block*, int> Block2Pos;

class MutationsData : public ThreadData {
public:
    Genome2Str genome2str;
    Block2Pos block2start;
    Block2Pos block2stop;
};

ThreadData* MutationsSequences::before_thread_impl() const {
    return new MutationsData;
}

typedef std::set<int> Positions;

static void add_positions(Positions& positions,
                          const Mutation& m,
                          int distance, int block_length) {
    int start = std::max(m.start - distance, 0);
    int stop = std::min(m.stop + distance, block_length - 1);
    for (int pos = start; pos <= stop; pos++) {
        positions.insert(pos);
    }
}

bool MutationsSequences::process_block_impl(Block* block,
        ThreadData* data) const {
    Positions positions;
    int distance = opt_value("mutation-distance").as<int>();
    int block_length = block->alignment_length();
    print_mutations_->find_mutations(block,
                                     boost::bind(add_positions,
                                             boost::ref(positions),
                                             _1, distance,
                                             block_length));
    MutationsData* d;
    d = boost::polymorphic_downcast<MutationsData*>(data);
    Genome2Str& genome2str = d->genome2str;
    Block2Pos& block2start = d->block2start;
    Block2Pos& block2stop = d->block2stop;
    if (genome2str.empty()) {
        block2start[block] = 0;
    } else {
        // length of one of sequences
        block2start[block] = genome2str.begin()->second.size();
    }
    int new_block_length = positions.size();
    block2stop[block] = block2start[block] + new_block_length - 1;
    BOOST_FOREACH (Fragment* f, *block) {
        std::string genome = f->seq()->genome();
        std::string& s = genome2str[genome];
        BOOST_ASSERT_MSG(s.size() == block2start[block],
                         "Forgot Steam --exact=1?");
        BOOST_FOREACH (int pos, positions) {
            // set is ordered
            char c = f->alignment_at(pos);
            if (c == '\0') {
                c = 'N';
            }
            s += c;
        }
    }
}

bool MutationsSequences::after_thread_impl(ThreadData* data) const {
    MutationsData* d;
    d = boost::polymorphic_downcast<MutationsData*>(data);
    BlockSet& bs = *block_set();
    if (bs.seqs().empty()) {
        BOOST_FOREACH (const Genome2Str::value_type& g_and_str,
                      d->genome2str) {
            const std::string& genome = g_and_str.first;
            SequencePtr seq = create_sequence(this);
            seq->set_name(genome);
            bs.add_sequence(seq);
        }
    }
    int shift = 0;
    BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
        const std::string& genome = seq->name();
        const std::string& str = d->genome2str[genome];
        shift = seq->size();
        seq->push_back(str);
    }
    BOOST_FOREACH (const Block2Pos::value_type& b_and_start,
                  d->block2start) {
        Block* orig_block = b_and_start.first;
        int start = b_and_start.second;
        int stop = d->block2stop[orig_block];
        int fragment_start = shift + start;
        int fragment_stop = shift + stop;
        Block* new_block = new Block;
        bs.insert(new_block);
        new_block->set_name(orig_block->name());
        BOOST_FOREACH (SequencePtr seq, bs.seqs()) {
            new_block->insert(new Fragment(seq, fragment_start,
                                           fragment_stop));
        }
    }
}

const char* MutationsSequences::name_impl() const {
    return "Create a seq. per genome of mutations and map blocks";
}

}

