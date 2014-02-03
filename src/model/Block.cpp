/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <cctype>
#include <map>
#include <set>
#include <algorithm>
#include "boost-xtime.hpp"
#include <boost/foreach.hpp>
#include <boost/pool/singleton_pool.hpp>
#include <boost/lexical_cast.hpp>

#include "Block.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "block_stat.hpp"
#include "block_hash.hpp"
#include "rand_name.hpp"
#include "throw_assert.hpp"

namespace bloomrepeats {

const int BLOCK_RAND_NAME_SIZE = 8;

Block::Block():
    name_(BLOCK_RAND_NAME_SIZE, '0'),
    weak_(false)
{ }

Block::Block(const std::string& name):
    weak_(false) {
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

void Block::insert(Fragment* fragment) {
    fragments_.push_back(fragment);
    if (!weak() || !fragment->block_raw_ptr()) {
        fragment->set_block(this);
    }
}

void Block::erase(Fragment* fragment) {
    Impl::iterator it = std::find(begin(), end(), fragment);
    BOOST_ASSERT(it != end());
    fragments_.erase(it);
    if (fragment->block_raw_ptr() == this) {
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

bool Block::has(Fragment* fragment) const {
    return std::find(begin(), end(), fragment) != end();
}

void Block::clear() {
    BOOST_FOREACH (Fragment* fragment, *this) {
        if (fragment->block_raw_ptr() == this) {
            fragment->set_block(0);
            delete fragment;
        }
    }
    fragments_.clear();
}

void Block::swap(Block& other) {
    fragments_.swap(other.fragments_);
    name_.swap(other.name_);
    std::swap(weak_, other.weak_);
    if (!this->weak()) {
        BOOST_FOREACH (Fragment* f, *this) {
            f->set_block(this);
        }
    }
    if (!other.weak()) {
        BOOST_FOREACH (Fragment* f, other) {
            f->set_block(&other);
        }
    }
}

Fragment* Block::front() const {
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
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        return f1->length() < f2->length();
    }
} fcl;

size_t Block::alignment_length() const {
    size_t result = 0;
    BOOST_FOREACH (Fragment* f, *this) {
        result = std::max(result, f->alignment_length());
    }
    return result;
}

float Block::identity() const {
    AlignmentStat al_stat;
    make_stat(al_stat, this);
    return block_identity(al_stat);
}

char Block::consensus_char(int pos, char gap) const {
    enum {
        A, T, C, G, N
    };
    const int LETTERS_NUMBER = 5;
    int freq[LETTERS_NUMBER] = {0, 0, 0, 0, 0};
    BOOST_FOREACH (Fragment* f, *this) {
        char c = f->alignment_at(pos);
        if (c == 'A') {
            freq[A] += 1;
        } else if (c == 'T') {
            freq[T] += 1;
        } else if (c == 'G') {
            freq[G] += 1;
        } else if (c == 'C') {
            freq[C] += 1;
        } else if (c == 'N') {
            freq[N] += 1;
        }
    }
    int max_freq = 0;
    for (int letter = 0; letter < LETTERS_NUMBER; letter++) {
        if (freq[letter] > max_freq) {
            max_freq = freq[letter];
        }
    }
    if (freq[A] == max_freq) {
        return 'A';
    } else if (freq[T] == max_freq) {
        return 'T';
    } else if (freq[G] == max_freq) {
        return 'G';
    } else if (freq[C] == max_freq) {
        return 'C';
    } else if (freq[N] == max_freq) {
        return 'N';
    } else {
        return gap;
    }
}

void Block::consensus(std::ostream& o, char gap) const {
    if (!empty() && !front()->row()) {
        Fragment* longest = front();
        BOOST_FOREACH (Fragment* f, *this) {
            BOOST_ASSERT_MSG(!f->row(), "Alignment rows are set to some of "
                             "fragments of block, being not set for other");
            if (f->length() > longest->length()) {
                longest = f;
            }
        }
        longest->print_contents(o);
    } else {
        int length = alignment_length();
        for (size_t pos = 0; pos < length; pos++) {
            o << consensus_char(pos, gap);
        }
    }
}

std::string Block::consensus_string(char gap) const {
    std::stringstream ss;
    consensus(ss, gap);
    return ss.str();
}

int Block::match(Block* one, Block* another) {
    if (one->size() != another->size()) {
        return 0;
    }
    bool all_match = true;
    bool all_match_inversed = true;
    typedef std::map<int, int> OriCount;
    typedef std::map<Sequence*, OriCount> Seq2Ori;
    Seq2Ori seq2ori, seq2ori_other;
    BOOST_FOREACH (Fragment* fragment, *one) {
        seq2ori[fragment->seq()][fragment->ori()] += 1;
    }
    BOOST_FOREACH (Fragment* fragment, *another) {
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

void Block::inverse() {
    BOOST_FOREACH (Fragment* fragment, *this) {
        fragment->inverse();
    }
}

void Block::patch(const FragmentDiff& diff) {
    BOOST_FOREACH (Fragment* fragment, *this) {
        fragment->patch(diff);
    }
}

Block* Block::split(size_t new_length) {
    Block* result = new Block();
    BOOST_FOREACH (Fragment* fragment, *this) {
        Fragment* new_fragment = fragment->split(new_length);
        if (new_fragment) {
            result->insert(new_fragment);
        }
    }
    return result;
}

void Block::find_place() {
    BOOST_FOREACH (Fragment* fragment, *this) {
        fragment->find_place();
    }
}

int Block::max_shift_end(int max_overlap) const {
    int result = INT_MAX;
    BOOST_FOREACH (Fragment* f, *this) {
        result = std::min(result, f->max_shift_end(max_overlap));
    }
    return result;
}

size_t Block::common_positions(const Fragment& fragment) const {
    size_t result = 0;
    BOOST_FOREACH (Fragment* f, *this) {
        result += f->common_positions(fragment);
    }
    return result;
}

void Block::merge(Block* other) {
    BOOST_ASSERT(!weak());
    BOOST_ASSERT(!other->weak());
    typedef std::map<Fragment, Fragment*> F2F;
    F2F f2f;
    std::vector<Fragment*> this_copy(begin(), end());
    BOOST_FOREACH (Fragment* f, this_copy) {
        if (f2f.find(*f) != f2f.end()) {
            delete f;
        } else {
            f2f[*f] = f;
        }
    }
    fragments_.clear();
    bool inverse_needed = false;
    BOOST_FOREACH (Fragment* f, *other) {
        f->inverse();
        if (f2f.find(*f) != f2f.end()) {
            inverse_needed = true;
        }
        f->inverse();
    }
    if (inverse_needed) {
        other->inverse();
    }
    std::vector<Fragment*> other_copy(other->begin(), other->end());
    BOOST_FOREACH (Fragment* f, other_copy) {
        if (f2f.find(*f) == f2f.end()) {
            f2f[*f] = f;
        } else {
            delete f;
        }
    }
    other->fragments_.clear();
    BOOST_FOREACH (F2F::value_type& f_and_ptr, f2f) {
        Fragment* f = f_and_ptr.second;
        insert(f);
    }
}

void Block::set_name(const std::string& name) {
#ifndef NDEBUG
    BOOST_FOREACH (char c, name) {
        BOOST_ASSERT(isgraph(c));
    }
#endif
    name_ = name;
}

void Block::set_random_name() {
    set_name(rand_name(BLOCK_RAND_NAME_SIZE));
}

void Block::set_name_from_fragments() {
    const char* const NAME_ABC = "0123456789abcdef";
    const int NAME_ABC_SIZE = 16;
    uint32_t a = block_hash(this);
    name_.resize(8);
    for (int byte_index = 0; byte_index < 4; byte_index++) {
        int byte = 0xFF & (a >> (8 * (3 - byte_index)));
        name_[byte_index * 2] = NAME_ABC[byte >> 4];
        name_[byte_index * 2 + 1] = NAME_ABC[byte & 0x0F];
    }
}

void Block::set_weak(bool weak) {
    if (this->weak() && !weak) {
        BOOST_FOREACH (Fragment* fragment, *this) {
            BOOST_ASSERT(fragment->block());
            if (fragment->block() != this) {
                if (!fragment->block()->weak()) {
                    fragment->block()->set_weak(true);
                }
                fragment->set_block(this);
            }
        }
    }
    weak_ = weak;
}

static struct FragmentCompareId {
    bool operator()(const Fragment* f1, const Fragment* f2) const {
        return f1->id() < f2->id();
    }
} fci;

std::ostream& operator<<(std::ostream& o, const Block& b) {
    std::vector<Fragment*> fragments(b.begin(), b.end());
    std::sort(fragments.begin(), fragments.end(), fci);
    BOOST_FOREACH (Fragment* f, fragments) {
        o << '>';
        f->print_header(o, &b);
        o << '\n';
        f->print_contents(o);
        o << '\n';
    }
    return o;
}

}

