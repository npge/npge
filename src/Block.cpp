/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <climits>
#include <cctype>
#include <map>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/pool/singleton_pool.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread/mutex.hpp>

#include "Block.hpp"
#include "Fragment.hpp"
#include "PairAligner.hpp"
#include "JoinApprover.hpp"

namespace bloomrepeats {

BlockPtr Block::create_new() {
    return new Block;
}

int block_number = 0;
boost::mutex next_block_mutex;

static int next_block_name() {
    boost::mutex::scoped_lock lock(next_block_mutex);
    block_number += 1;
    return block_number;
}

Block::Block():
    name_(boost::lexical_cast<std::string>(next_block_name()))
{ }

Block::Block(const std::string& name) {
    set_name(name);
}

Block::~Block() {
    clear();
}

class BlockTag;

typedef boost::singleton_pool<BlockTag, sizeof(Block)> BlockPool;

void* Block::operator new(size_t /* x */) {
    return BlockPool::malloc();
}

void Block::operator delete(void* ptr) {
    BlockPool::free(ptr);
}

BlockPtr Block::clone() const {
    BlockPtr result = create_new();
    BOOST_FOREACH (FragmentPtr f, *this) {
        result->insert(Fragment::create_new(*f));
    }
    return result;
}

void Block::insert(FragmentPtr fragment) {
#ifndef NDEBUG
    BOOST_FOREACH (FragmentPtr f, *this) {
        BOOST_ASSERT(*fragment != *f);
    }
#endif
    fragments_.push_back(fragment);
    fragment->set_block(this);
}

void Block::erase(FragmentPtr fragment) {
    Impl::iterator it = std::find(begin(), end(), fragment);
    BOOST_ASSERT(it != end());
    fragments_.erase(it);
    if (fragment->block_raw_ptr()) {
        fragment->set_block(0);
        delete fragment;
    }
}

size_t Block::size() const {
    return fragments_.size();
}

bool Block::empty() const {
    return fragments_.empty();
}

bool Block::has(FragmentPtr fragment) const {
    return std::find(begin(), end(), fragment) != end();
}

void Block::clear() {
    BOOST_FOREACH (FragmentPtr fragment, *this) {
        if (fragment->block_raw_ptr()) {
            fragment->set_block(0);
            delete fragment;
        }
    }
    fragments_.clear();
}

FragmentPtr Block::front() const {
    return empty() ? 0 : *(begin());
}

Block::Impl::iterator Block::begin() {
    return fragments_.begin();
}

Block::Impl::const_iterator Block::begin() const {
    return fragments_.begin();
}

Block::Impl::iterator Block::end() {
    return fragments_.end();
}

Block::Impl::const_iterator Block::end() const {
    return fragments_.end();
}

static struct FragmentCompareLength {
    bool operator()(const FragmentPtr& f1, const FragmentPtr& f2) const {
        return f1->length() < f2->length();
    }
} fcl;

float Block::identity() const {
    size_t total = (*std::max_element(begin(), end(), fcl))->length();
    size_t equal = 0;
    size_t min_length = (*std::min_element(begin(), end(), fcl))->length();
    for (size_t pos = 0; pos < min_length; pos++) {
        char c = 0;
        BOOST_FOREACH (FragmentPtr f, *this) {
            if (c == 0) {
                c = f->at(pos);
            } else if (c != f->at(pos)) {
                c = -1;
                break;
            }
        }
        if (c != -1) {
            equal += 1;
        }
    }
    return float(equal) / float(total);
}

int Block::match(BlockPtr one, BlockPtr another) {
    if (one->size() != another->size()) {
        return 0;
    }
    bool all_match = true;
    bool all_match_inversed = true;
    typedef std::map<int, int> OriCount;
    typedef std::map<Sequence*, OriCount> Seq2Ori;
    Seq2Ori seq2ori, seq2ori_other;
    BOOST_FOREACH (FragmentPtr fragment, *one) {
        seq2ori[fragment->seq()][fragment->ori()] += 1;
    }
    BOOST_FOREACH (FragmentPtr fragment, *another) {
        seq2ori_other[fragment->seq()][fragment->ori()] += 1;
    }
    BOOST_FOREACH (Seq2Ori::value_type& seq_and_ori, seq2ori) {
        Sequence* seq = seq_and_ori.first;
        OriCount& ori_count = seq_and_ori.second;
        Seq2Ori::iterator it = seq2ori_other.find(seq);
        if (it == seq2ori_other.end()) {
            return 0;
        }
        OriCount& ori_count_other = it->second;
        for (int ori = -1; ori <= 1; ori += 2) {
            if (ori_count[ori] != ori_count_other[ori]) {
                all_match = false;
            }
            if (ori_count[ori] != ori_count_other[-ori]) {
                all_match_inversed = false;
            }
        }
        if (!all_match && !all_match_inversed) {
            return 0;
        }
    }
    BOOST_ASSERT(all_match || all_match_inversed);
    return all_match ? 1 : -1;
}

void Block::filter(int min_fragment_length) {
    std::vector<FragmentPtr> block_copy(begin(), end());
    BOOST_FOREACH (FragmentPtr fragment, block_copy) {
        if (!fragment->valid() || fragment->length() < min_fragment_length) {
            fragment->disconnect();
            erase(fragment);
        }
    }
}

int Block::can_join(BlockPtr one, BlockPtr another) {
    bool all[3] = {true, false, true};
    for (int ori = 1; ori >= -1; ori -= 2) {
        BOOST_FOREACH (FragmentPtr f, *one) {
            FragmentPtr f1 = f->logical_neighbor(ori);
            if (!f1 || f1->block() != another || !Fragment::can_join(f, f1)) {
                all[ori + 1] = false;
                break;
            }
        }
        if (all[ori + 1]) {
            break;
        }
    }
    int result = all[1 + 1] ? 1 : all[-1 + 1] ? -1 : 0;
    BOOST_ASSERT(!(result && !match(one, another)));
    return result;
}

BlockPtr Block::join(BlockPtr one, BlockPtr another, int logical_ori) {
    BOOST_ASSERT(can_join(one, another) == logical_ori);
    BlockPtr result = create_new();
    std::set<FragmentPtr> to_delete;
    BOOST_FOREACH (FragmentPtr f, *one) {
        FragmentPtr f1 = f->logical_neighbor(logical_ori);
        BOOST_ASSERT(f1);
        result->insert(Fragment::join(f, f1));
        to_delete.insert(f);
        to_delete.insert(f1);
    }
    BOOST_FOREACH (FragmentPtr f, to_delete) {
        delete f;
    }
    return result;
}

BlockPtr Block::try_join(BlockPtr one, BlockPtr another, JoinApprover* ja) {
    BlockPtr result = 0;
    int match_ori = match(one, another);
    if (match_ori == -1) {
        another->inverse();
    }
    if (match_ori) {
        int logical_ori = can_join(one, another);
        if (logical_ori && (!ja || ja->can_join_blocks(one, another))) {
            result = join(one, another, logical_ori);
        }
    }
    return result;
}

void Block::inverse() {
    BOOST_FOREACH (FragmentPtr fragment, *this) {
        fragment->inverse();
    }
}

void Block::patch(const FragmentDiff& diff) {
    BOOST_FOREACH (FragmentPtr fragment, *this) {
        fragment->patch(diff);
    }
}

BlockPtr Block::split(size_t new_length) {
    BlockPtr result = Block::create_new();
    BOOST_FOREACH (FragmentPtr fragment, *this) {
        FragmentPtr new_fragment = fragment->split(new_length);
        if (new_fragment) {
            result->insert(new_fragment);
        }
    }
    return result;
}

void Block::find_place() {
    BOOST_FOREACH (FragmentPtr fragment, *this) {
        fragment->find_place();
    }
}

int Block::max_shift_end(int max_overlap) const {
    int result = INT_MAX;
    BOOST_FOREACH (FragmentPtr f, *this) {
        result = std::min(result, f->max_shift_end(max_overlap));
    }
    return result;
}

void Block::expand(PairAligner* aligner, int batch, int ori, int max_overlap) {
    aligner = aligner ? : PairAligner::default_aligner();
    if (ori == 1) {
        if (size() >= 2) {
            expand_end(*aligner, batch, max_overlap);
        }
    } else if (ori == -1) {
        inverse();
        expand(aligner, batch, /* ori */ 1, max_overlap);
        inverse();
    } else { /* ori = 0 */
        expand(aligner, batch, /* ori */ 1, max_overlap);
        expand(aligner, batch, /* ori */ -1, max_overlap);
    }
}

size_t Block::common_positions(const Fragment& fragment) {
    size_t result = 0;
    BOOST_FOREACH (FragmentPtr f, *this) {
        result += f->common_positions(fragment);
    }
    return result;
}

bool Block::expand_by_fragments(PairAligner* aligner, int batch) {
    bool result = false;
    aligner = aligner ? : PairAligner::default_aligner();
    std::set<BlockPtr> visited;
    BOOST_FOREACH (FragmentPtr f, std::vector<FragmentPtr>(begin(), end())) {
        for (int ori = 1; ori >= -1; ori -= 2) {
            FragmentPtr neighbor = f->neighbor(ori);
            if (neighbor) {
                FragmentDiff diff = neighbor->diff_to(*f);
                BlockPtr block = neighbor->block();
                if (block && block != this &&
                        visited.find(block) == visited.end()) {
                    visited.insert(block);
                    BOOST_FOREACH (FragmentPtr fn, *block) {
                        Fragment candidate;
                        candidate.apply_coords(*fn);
                        candidate.patch(diff);
                        if (candidate.valid() && !common_positions(candidate) &&
                                f->aligned(candidate, aligner, batch)) {
                            FragmentPtr new_f = Fragment::create_new();
                            new_f->apply_coords(candidate);
                            insert(new_f);
                            new_f->find_place(fn);
                            result = true;
                        }
                    }
                }
            }
        }
    }
    return result;
}

void Block::merge(BlockPtr other) {
    typedef std::map<Fragment, FragmentPtr> F2F;
    F2F f2f;
    BOOST_FOREACH (FragmentPtr f, *this) {
        BOOST_ASSERT(f2f.find(*f) == f2f.end());
        f2f[*f] = f;
    }
    fragments_.clear();
    bool inverse_needed = false;
    BOOST_FOREACH (FragmentPtr f, *other) {
        f->inverse();
        if (f2f.find(*f) != f2f.end()) {
            inverse_needed = true;
        }
        f->inverse();
    }
    if (inverse_needed) {
        other->inverse();
    }
    std::vector<FragmentPtr> other_copy(other->begin(), other->end());
    BOOST_FOREACH (FragmentPtr f, other_copy) {
        if (f2f.find(*f) == f2f.end()) {
            f2f[*f] = f;
        } else {
            delete f;
        }
    }
    other->fragments_.clear();
    BOOST_FOREACH (F2F::value_type& f_and_ptr, f2f) {
        FragmentPtr f = f_and_ptr.second;
        insert(f);
    }
}

void Block::set_name(const std::string& name) {
    BOOST_ASSERT(name.length() >= 1 && name.length() <= 40);
#ifndef NDEBUG
    BOOST_FOREACH (char c, name) {
        BOOST_ASSERT(isalnum(c));
    }
#endif
    name_ = name;
}

void Block::expand_end(PairAligner& aligner, int batch, int max_overlap) {
    std::vector<int> main_end(size() - 1), o_end(size() - 1);
    FragmentPtr main_f = fragments_.back();
    while (true) {
        int max_shift = max_shift_end(max_overlap);
        if (max_shift <= 0) {
            break;
        }
        int shift = std::min(batch, max_shift);
        std::string main_str = main_f->substr(-1, main_f->length() - 1 + shift);
        aligner.set_first(main_str.c_str(), main_str.size());
        for (int i = 0; i < fragments_.size() - 1; i++) {
            FragmentPtr o_f = fragments_[i];
            std::string o_str = o_f->substr(-1, o_f->length() - 1 + shift);
            aligner.set_second(o_str.c_str(), o_str.size());
            aligner.align(main_end[i], o_end[i]);
        }
        int min_end = *std::min_element(main_end.begin(), main_end.end());
        main_f->shift_end(min_end);
        for (int i = 0; i < fragments_.size() - 1; i++) {
            FragmentPtr o_f = fragments_[i];
            int delta = main_end[i] - min_end;
            o_f->shift_end(o_end[i] - delta);
        }
        const float MIN_ACCEPTED = 0.5;
        if (min_end < batch * MIN_ACCEPTED) {
            break;
        }
    }
}

std::ostream& operator<<(std::ostream& o, const Block& b) {
    BOOST_FOREACH (FragmentPtr f, b) {
        o << *f << std::endl;
    }
    return o;
}

}

