/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
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
#include "Exception.hpp"
#include "throw_assert.hpp"
#include "block_hash.hpp"
#include "global.hpp"

namespace npge {

typedef std::map<std::string, SequencePtr> Name2Seq;
typedef std::map<std::string, BSA> Name2BSA;

struct BlockSet::I {
    BlockSet::Impl blocks_;
    std::set<SequencePtr> seqs_;
    Name2BSA bsas_;
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

std::string BlockSet::block_from_description(const std::string& description) {
    return extract_value(description, "block");
}

void BlockSet::insert(Block* block) {
    Impl& blocks = impl_->blocks_;
    ASSERT_TRUE(blocks.find(block) == blocks.end());
    blocks.insert(block);
}

void BlockSet::erase(Block* block) {
    detach(block);
    delete block;
}

void BlockSet::detach(Block* block) {
    impl_->blocks_.erase(block);
}

int BlockSet::size() const {
    return impl_->blocks_.size();
}

bool BlockSet::empty() const {
    return impl_->blocks_.empty();
}

bool BlockSet::has(const Block* block) const {
    Block* b = const_cast<Block*>(block);
    return impl_->blocks_.find(b) != impl_->blocks_.end();
}

Block* BlockSet::find_block(const std::string& name) const {
    BOOST_FOREACH (Block* b, *this) {
        if (b->name() == name) {
            return b;
        }
    }
    return 0;
}

void BlockSet::clear() {
    clear_blocks();
    clear_seqs();
    clear_bsas();
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

BSA& BlockSet::bsa(const std::string& bsa_name) {
    return impl_->bsas_[bsa_name];
}

const BSA& BlockSet::bsa(const std::string& bsa_name) const {
    Name2BSA::const_iterator it = impl_->bsas_.find(bsa_name);
    if (it == impl_->bsas_.end()) {
        throw Exception("No such bsa: " + bsa_name);
    } else {
        return it->second;
    }
}

Strings BlockSet::bsas() const {
    Strings bsa_names;
    BOOST_FOREACH (const Name2BSA::value_type& n_b, impl_->bsas_) {
        bsa_names.push_back(n_b.first);
    }
    return bsa_names;
}

bool BlockSet::has_bsa(const std::string& bsa_name) const {
    return impl_->bsas_.find(bsa_name) != impl_->bsas_.end();
}

void BlockSet::remove_bsa(const std::string& bsa_name) {
    impl_->bsas_.erase(bsa_name);
}

void BlockSet::clear_bsas() {
    impl_->bsas_.clear();
}

void BlockSet::swap(BlockSet& other) {
    impl_->blocks_.swap(other.impl_->blocks_);
    impl_->seqs_.swap(other.impl_->seqs_);
    impl_->bsas_.swap(other.impl_->bsas_);
}

BlockSetPtr BlockSet::clone() const {
    BlockSetPtr result = new_bs();
    copy(*result);
    return result;
}

void BlockSet::copy(BlockSet& target) const {
    BOOST_FOREACH (Block* block, *this) {
        target.insert(block->clone());
    }
    target.add_sequences(seqs());
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

bool BlockSet::operator==(const BlockSet& other) const {
    return blockset_hash(*this) == blockset_hash(other);
}

std::istream& operator>>(std::istream& input, BlockSet& block_set) {
    BlockSetFastaReader reader(block_set, input,
                               COMPACT_ROW,
                               COMPACT_LOW_N_SEQUENCE);
    reader.run();
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

