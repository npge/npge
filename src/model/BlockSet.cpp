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
#include "boost-xtime.hpp"
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/thread/shared_mutex.hpp>

#include "BlockSet.hpp"
#include "read_block_set.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "AlignmentRow.hpp"
#include "key_value.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

typedef std::map<std::string, SequencePtr> Name2Seq;

struct BlockSet::I {
    BlockSet::Impl blocks_;
    std::set<SequencePtr> seqs_;
    Name2Seq name2seq_;
    boost::shared_mutex name2seq_mutex_;
};

BlockSet::BlockSet() {
    impl_ = new I;
}

BlockSet::~BlockSet() {
    clear();
    delete impl_;
    impl_ = 0;
}

void BlockSet::add_sequence(SequencePtr seq) {
    impl_->seqs_.insert(seq);
}

void BlockSet::add_sequences(const std::vector<SequencePtr>& sequences) {
    BOOST_FOREACH (const SequencePtr& seq, sequences) {
        add_sequence(seq);
    }
}

std::vector<SequencePtr> BlockSet::seqs() const {
    return std::vector<SequencePtr>(impl_->seqs_.begin(), impl_->seqs_.end());
}

void BlockSet::remove_sequence(const SequencePtr& seq) {
    impl_->seqs_.erase(seq);
}

SequencePtr BlockSet::seq_from_name(const std::string& name) const {
    boost::upgrade_lock<boost::shared_mutex> lock(impl_->name2seq_mutex_);
    Name2Seq::const_iterator it = impl_->name2seq_.find(name);
    if (it == impl_->name2seq_.end()) {
        boost::upgrade_to_unique_lock<boost::shared_mutex> unique_lock(lock);
        BOOST_FOREACH (SequencePtr seq, impl_->seqs_) {
            impl_->name2seq_[seq->name()] = seq;
        }
        it = impl_->name2seq_.find(name);
    }
    if (it != impl_->name2seq_.end()) {
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
    return extract_value(description, "block");
}

void BlockSet::insert(Block* block) {
#ifndef NDEBUG
    BOOST_FOREACH (Block* b, *this) {
        BOOST_ASSERT(block != b);
    }
#endif
    impl_->blocks_.insert(block);
}

void BlockSet::erase(Block* block) {
    detach(block);
    delete block;
}

void BlockSet::detach(Block* block) {
    impl_->blocks_.erase(block);
}

size_t BlockSet::size() const {
    return impl_->blocks_.size();
}

bool BlockSet::empty() const {
    return impl_->blocks_.empty();
}

bool BlockSet::has(const Block* block) const {
    Block* b = const_cast<Block*>(block);
    return impl_->blocks_.find(b) != impl_->blocks_.end();
}

void BlockSet::clear() {
    clear_blocks();
    clear_seqs();
}

void BlockSet::clear_blocks() {
    BOOST_FOREACH (Block* block, *this) {
        delete block;
    }
    impl_->blocks_.clear();
}

void BlockSet::clear_seqs() {
    impl_->seqs_.clear();
}

void BlockSet::swap(BlockSet& other) {
    impl_->blocks_.swap(other.impl_->blocks_);
    impl_->seqs_.swap(other.impl_->seqs_);
    impl_->name2seq_.swap(other.impl_->name2seq_);
}

Block* BlockSet::front() const {
    return empty() ? 0 : *(begin());
}

BlockSet::iterator BlockSet::begin() {
    return impl_->blocks_.begin();
}

BlockSet::const_iterator BlockSet::begin() const {
    return impl_->blocks_.begin();
}

BlockSet::iterator BlockSet::end() {
    return impl_->blocks_.end();
}

BlockSet::const_iterator BlockSet::end() const {
    return impl_->blocks_.end();
}

std::istream& operator>>(std::istream& input, BlockSet& block_set) {
    BlockSetFastaReader reader(block_set, input,
                               COMPACT_ROW, COMPACT_SEQUENCE);
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

