/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "BlockSet.hpp"
#include "read_block_set.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "AlignmentRow.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

BlockSet::~BlockSet() {
    clear();
}

void BlockSet::add_sequence(SequencePtr seq) {
    seqs_.insert(seq);
}

void BlockSet::add_sequences(const std::vector<SequencePtr>& sequences) {
    BOOST_FOREACH (const SequencePtr& seq, sequences) {
        add_sequence(seq);
    }
}

std::vector<SequencePtr> BlockSet::seqs() const {
    return std::vector<SequencePtr>(seqs_.begin(), seqs_.end());
}

SequencePtr BlockSet::seq_from_name(const std::string& name) const {
    boost::upgrade_lock<boost::shared_mutex> lock(name2seq_mutex_);
    Name2Seq::const_iterator it = name2seq_.find(name);
    if (it == name2seq_.end()) {
        boost::upgrade_to_unique_lock<boost::shared_mutex> unique_lock(lock);
        BOOST_FOREACH (SequencePtr seq, seqs_) {
            name2seq_[seq->name()] = seq;
        }
        it = name2seq_.find(name);
    }
    if (it != name2seq_.end()) {
        return it->second;
    } else {
        return SequencePtr();
    }
}

Fragment* BlockSet::fragment_from_id(const std::string& id) const {
    if (id.empty()) {
        return 0;
    }
    size_t u1 = id.find('_');
    if (u1 == std::string::npos) {
        return 0;
    }
    std::string seq_name = id.substr(0, u1);
    if (seq_name.empty()) {
        return 0;
    }
    SequencePtr seq = seq_from_name(seq_name);
    if (!seq) {
        return 0;
    }
    size_t u2 = id.find('_', u1 + 1);
    if (u2 == std::string::npos) {
        return 0;
    }
    std::string begin_pos_str = id.substr(u1 + 1, u2 - u1 - 1);
    size_t begin_pos = boost::lexical_cast<size_t>(begin_pos_str);
    std::string last_pos_str = id.substr(u2 + 1);
    size_t last_pos = boost::lexical_cast<size_t>(last_pos_str);
    Fragment* f = new Fragment(seq);
    f->set_ori(begin_pos <= last_pos ? 1 : -1);
    f->set_begin_pos(begin_pos);
    f->set_last_pos(last_pos);
    return f;
}

std::string BlockSet::block_from_description(const std::string& description) {
    size_t block_pos = description.find("block=");
    if (block_pos == std::string::npos) {
        return "";
    }
    size_t block_name_start = block_pos + std::string("block=").size();
    size_t space_pos = description.find(' ', block_name_start); // or npos
    return description.substr(block_name_start, space_pos - block_name_start);
}

void BlockSet::insert(Block* block) {
#ifndef NDEBUG
    BOOST_FOREACH (Block* b, *this) {
        BOOST_ASSERT(block != b);
    }
#endif
    blocks_.insert(block);
}

void BlockSet::erase(Block* block) {
    blocks_.erase(block);
    delete block;
}

size_t BlockSet::size() const {
    return blocks_.size();
}

bool BlockSet::empty() const {
    return blocks_.empty();
}

bool BlockSet::has(const Block* block) const {
    return blocks_.find(const_cast<Block*>(block)) != blocks_.end();
}

void BlockSet::clear() {
    BOOST_FOREACH (Block* block, *this) {
        delete block;
    }
    blocks_.clear();
}

void BlockSet::swap(BlockSet& other) {
    blocks_.swap(other.blocks_);
    seqs_.swap(other.seqs_);
    name2seq_.swap(other.name2seq_);
}

Block* BlockSet::front() const {
    return empty() ? 0 : *(begin());
}

BlockSet::iterator BlockSet::begin() {
    return blocks_.begin();
}

BlockSet::const_iterator BlockSet::begin() const {
    return blocks_.begin();
}

BlockSet::iterator BlockSet::end() {
    return blocks_.end();
}

BlockSet::const_iterator BlockSet::end() const {
    return blocks_.end();
}

std::istream& operator>>(std::istream& input, BlockSet& block_set) {
    BlockSetFastaReader reader(block_set, input,
                               /* keep_alignment */ true, COMPACT_ROW);
    reader.read_all_sequences();
    return input;
}

static struct BlockCompareName {
    bool operator()(const Block* b1, const Block* b2) const {
        typedef boost::tuple<int, const std::string&> Tie;
        return Tie(-b1->size(), b1->name()) < Tie(-b2->size(), b2->name());
    }
} bcn;

std::ostream& operator<<(std::ostream& o, const BlockSet& block_set) {
    std::vector<Block*> blocks(block_set.begin(), block_set.end());
    std::sort(blocks.begin(), blocks.end(), bcn);
    BOOST_FOREACH (Block* block, blocks) {
        o << *block;
        o << std::endl; // empty line
    }
    return o;
}

BlockSetPtr new_bs() {
    return boost::make_shared<BlockSet>();
}

}

