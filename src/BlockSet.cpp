/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <vector>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "BlockSet.hpp"
#include "Block.hpp"
#include "Fragment.hpp"
#include "Sequence.hpp"
#include "PairAligner.hpp"

namespace bloomrepeats {

void BlockSet::add_sequence(SequencePtr seq) {
    seqs_.push_back(seq);
}

void BlockSet::insert(BlockPtr block) {
#ifndef NDEBUG
    BOOST_FOREACH (BlockPtr b, *this) {
        BOOST_ASSERT(block != b);
    }
#endif
    blocks_.insert(block);
}

void BlockSet::erase(BlockPtr block) {
    blocks_.erase(block);
}

size_t BlockSet::size() const {
    return blocks_.size();
}

bool BlockSet::empty() const {
    return blocks_.empty();
}

bool BlockSet::has(BlockPtr block) const {
    return blocks_.find(block) != blocks_.end();
}

void BlockSet::clear() {
    blocks_.clear();
}

BlockPtr BlockSet::front() const {
    return empty() ? BlockPtr() : *(begin());
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

static struct FragmentCompare {
    bool operator()(const FragmentPtr& f1, const FragmentPtr& f2) const {
        return *f1 < *f2;
    }
} fragment_compare;

void BlockSet::connect_fragments() {
    typedef std::vector<FragmentPtr> Fs;
    typedef std::map<Sequence*, Fs> Seq2Fs;
    Seq2Fs seq2fs;
    BOOST_FOREACH (const BlockPtr& block, *this) {
        BOOST_FOREACH (FragmentPtr fragment, *block) {
            seq2fs[fragment->seq()].push_back(fragment);
        }
    }
    BOOST_FOREACH (Seq2Fs::value_type& seq_and_fs, seq2fs) {
        Fs& fs = seq_and_fs.second;
        std::sort(fs.begin(), fs.end(), fragment_compare);
        for (int i = 1; i < fs.size(); i++) {
            Fragment::connect(fs[i - 1], fs[i]);
        }
    }
}

void BlockSet::filter(int min_fragment_length, int min_block_size) {
    std::vector<BlockPtr> block_set_copy(begin(), end());
    BOOST_FOREACH (const BlockPtr& block, block_set_copy) {
        block->filter(min_fragment_length);
        if (block->size() < min_block_size) {
            erase(block);
        }
    }
}

static struct BlockGreater {
    bool operator()(const BlockPtr& b1, const BlockPtr& b2) const {
        return b1->size() > b2->size();
    }
} block_greater;

static BlockPtr neighbour_block(const BlockPtr& b, int ori) {
    BlockPtr result;
    FragmentPtr f = b->front();
    if (f) {
        FragmentPtr neighbour_f = ori == 1 ? f->next() : f->prev();
        if (neighbour_f) {
            result = neighbour_f->block();
        }
    }
    return result;
}

void BlockSet::join(size_t max_gap) {
    std::vector<BlockPtr> bs(begin(), end());
    std::sort(bs.begin(), bs.end(), block_greater);
    BOOST_FOREACH (BlockPtr block, bs) {
        if (has(block)) {
            for (int ori = -1; ori <= 1; ori += 2) {
                while (BlockPtr other_block = neighbour_block(block, ori)) {
                    BlockPtr new_block = Block::try_join(block, other_block,
                                                         max_gap);
                    if (new_block) {
                        erase(block);
                        erase(other_block);
                        insert(new_block);
                        block = new_block;
                    } else {
                        break;
                    }
                }
            }
        }
    }
}

void BlockSet::expand_blocks(PairAligner* aligner, int batch,
                             int ori, int max_overlap) {
    aligner = aligner ? : PairAligner::default_aligner();
    std::vector<BlockPtr> bs(begin(), end());
    std::sort(bs.begin(), bs.end(), block_greater);
    BOOST_FOREACH (BlockPtr block, bs) {
        block->expand(aligner, batch, ori, max_overlap);
    }
}

bool BlockSet::intersections() const {
    BOOST_FOREACH (BlockPtr block, *this) {
        BOOST_FOREACH (FragmentPtr fragment, *block) {
            for (int ori = -1; ori <= 1; ori += 2) {
                FragmentPtr neighbour = fragment->neighbour(ori);
                if (neighbour && fragment->common_positions(*neighbour)) {
                    return true;
                }
            }
        }
    }
    return false;
}

static struct BlockLess {
    bool operator()(const BlockPtr& b1, const BlockPtr& b2) const {
        return b1->size() < b2->size();
    }
} block_less;

typedef std::priority_queue<BlockPtr, std::vector<BlockPtr>, BlockLess> BQ;

static void treat_fragments(BlockSet* block_set, BQ& bs,
                            FragmentPtr x, FragmentPtr y) {
    BlockPtr x_block = x->block();
    BlockPtr y_block = y->block();
    if (x_block == y_block) {
        x_block->erase(x);
        return;
    }
    Fragment common = x->common_fragment(*y);
    BOOST_ASSERT(common.valid());
    if (*x == common && x->length() == y->length()) {
        BOOST_ASSERT(y_block);
        x->block()->merge(y_block);
        BOOST_ASSERT(y_block->empty());
        block_set->erase(y_block);
    } else {
        if (common == *x) {
            treat_fragments(block_set, bs, y, x);
        } else {
            size_t new_length;
            if (common.begin_pos() == x->begin_pos()) {
                new_length = common.length();
            } else {
                new_length = std::min(abs(x->begin_pos() - common.min_pos()),
                                      abs(x->begin_pos() - common.max_pos()));
            }
            BlockPtr new_block = x->block()->split(new_length);
            BOOST_ASSERT(new_block && !new_block->empty());
            bs.push(new_block);
            block_set->insert(new_block);
        }
    }
}

static bool treat_block(BlockSet* block_set, BQ& bs, const BlockPtr& block) {
    BOOST_FOREACH (FragmentPtr f, *block) {
        for (int ori = -1; ori <= 1; ori += 2) {
            FragmentPtr o_f = f->neighbour(ori);
            if (o_f && f->common_positions(*o_f)) {
                treat_fragments(block_set, bs, f, o_f);
                return true;
            }
        }
    }
    return false;
}

void BlockSet::resolve_intersections() {
    BQ bs(begin(), end(), block_less);
    while (!bs.empty()) {
        BlockPtr block = bs.top();
        bs.pop();
        while (has(block) && treat_block(this, bs, block))
        { }
    }
#ifndef NDEBUG
    BOOST_ASSERT(!intersections());
    connect_fragments();
    BOOST_ASSERT(!intersections());
#endif
}

bool BlockSet::expand_blocks_by_fragments(PairAligner* aligner, int batch) {
    aligner = aligner ? : PairAligner::default_aligner();
    bool result = false;
    BOOST_FOREACH (BlockPtr block, *this) {
        BOOST_ASSERT(block);
        result |= block->expand_by_fragments(aligner, batch);
    }
    return result;
}

static void try_new_block(BlockSet& set, const Fragment& f, int ori,
                          FragmentPtr* prev) {
    FragmentPtr n = f.neighbour(ori);
    FragmentPtr new_f = Fragment::create_new(f.seq());
    if (ori == -1) {
        new_f->set_min_pos(n ? n->max_pos() + 1 : 0);
        new_f->set_max_pos(f.min_pos() - 1);
    } else {
        new_f->set_min_pos(f.max_pos() + 1);
        new_f->set_max_pos(n ? n->min_pos() - 1 : f.seq()->size() - 1);
    }
    if (new_f->valid()) {
        if (*prev) {
            BOOST_ASSERT(!(*new_f < **prev));
            Fragment::connect(*prev, new_f);
        }
        *prev = new_f;
        BlockPtr block = Block::create_new();
        block->insert(new_f);
        set.insert(block);
    } else {
        delete new_f;
    }
}

BlockSetPtr BlockSet::rest() const {
    BlockSetPtr result = boost::make_shared<BlockSet>();
    result->seqs_ = seqs_;
    std::set<Sequence*> used;
    BOOST_FOREACH (const BlockPtr& block, *this) {
        BOOST_FOREACH (FragmentPtr f, *block) {
            Sequence* seq = f->seq();
            if (used.find(seq) == used.end()) {
                used.insert(seq);
                FragmentPtr prev = 0;
                while (FragmentPtr fr = f->neighbour(-1)) {
                    f = fr;
                }
                try_new_block(*result, *f, -1, &prev);
                while (FragmentPtr fr = f->neighbour(1)) {
                    f = fr;
                    try_new_block(*result, *f, -1, &prev);
                }
                try_new_block(*result, *f, 1, &prev);
            }
        }
    }
    return result;
}

void BlockSet::_read(std::istream& input,
                     const std::vector<SequencePtr>& seqs) {
    std::map<std::string, Sequence*> name2seq;
    BOOST_FOREACH (SequencePtr seq, seqs) {
        name2seq[seq->name()] = seq.get();
    }
    BlockPtr block = Block::create_new();
    for (std::string line; std::getline(input, line);) {
        boost::algorithm::trim(line);
        if (line.empty() && !block->empty()) {
            insert(block);
            block = Block::create_new();
        } else if (line.size() >= 1 && line[0] == '>') {
            size_t sp = line.find(' ');
            BOOST_ASSERT(line.size() >= 2);
            std::string name = line.substr(1, sp - 1);
            BOOST_ASSERT(!name.empty());
            size_t u1 = name.find('_');
            BOOST_ASSERT(u1 != std::string::npos);
            std::string seq_name = name.substr(0, u1);
            Sequence* seq = name2seq[seq_name];
            BOOST_ASSERT(seq);
            BOOST_ASSERT(!seq_name.empty());
            size_t u2 = name.find('_', u1 + 1);
            BOOST_ASSERT(u2 != std::string::npos);
            std::string begin_pos_str = name.substr(u1 + 1, u2 - u1 - 1);
            size_t begin_pos = boost::lexical_cast<size_t>(begin_pos_str);
            std::string last_pos_str = name.substr(u2 + 1);
            size_t last_pos = boost::lexical_cast<size_t>(last_pos_str);
            FragmentPtr f = Fragment::create_new(seq);
            f->set_ori(begin_pos < last_pos ? 1 : -1);
            f->set_begin_pos(begin_pos);
            f->set_last_pos(last_pos);
            block->insert(f);
        }
    }
    if (!block->empty()) {
        insert(block);
    }
}

std::ostream& operator<<(std::ostream& o, const BlockSet& block_set) {
    BOOST_FOREACH (BlockPtr block, block_set) {
        o << *block << std::endl;
    }
    return o;
}

}

