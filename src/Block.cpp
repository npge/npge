/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include <boost/thread/tss.hpp>

#include "Block.hpp"
#include "Fragment.hpp"
#include "PairAligner.hpp"

namespace bloomrepeats {

BlockPtr Block::create_new() {
    return boost::make_shared<Block>();
}

void Block::insert(FragmentPtr fragment) {
#ifndef NDEBUG
    BOOST_FOREACH (FragmentPtr f, *this) {
        BOOST_ASSERT(*fragment != *f);
    }
#endif
    fragments_.push_back(fragment);
    fragment->block_ = shared_from_this();
}

void Block::erase(FragmentPtr fragment) {
    Impl::iterator it = std::find(begin(), end(), fragment);
    BOOST_ASSERT(it != end());
    fragments_.erase(it);
    fragment->block_.reset();
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
        fragment->block_.reset();
    }
    fragments_.clear();
}

FragmentPtr Block::front() const {
    return empty() ? FragmentPtr() : *(begin());
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

int Block::match(const BlockPtr& one, const BlockPtr& another) {
    if (one->size() != another->size()) {
        return 0;
    }
    bool all_match = true;
    bool all_match_inversed = true;
    typedef std::map<int, int> OriCount;
    typedef std::map<SequencePtr, OriCount> Seq2Ori;
    Seq2Ori seq2ori, seq2ori_other;
    BOOST_FOREACH (const FragmentPtr& fragment, *one) {
        seq2ori[fragment->seq()][fragment->ori()] += 1;
    }
    BOOST_FOREACH (const FragmentPtr& fragment, *another) {
        seq2ori_other[fragment->seq()][fragment->ori()] += 1;
    }
    BOOST_FOREACH (Seq2Ori::value_type& seq_and_ori, seq2ori) {
        const SequencePtr& seq = seq_and_ori.first;
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
    BOOST_FOREACH (FragmentPtr fragment, *this) {
        fragment->inverse();
    }
}

boost::thread_specific_ptr<PairAligner> local_aligner_;

static PairAligner* local_aligner() {
    if (local_aligner_.get() == 0) {
        local_aligner_.reset(new PairAligner());
    }
    return local_aligner_.get();
}

void Block::expand(PairAligner* aligner, int batch, int ori) {
    aligner = aligner ? : local_aligner();
    if (ori == 1) {
        if (size() >= 2) {
            expand_end(*aligner, batch);
        }
    } else if (ori == -1) {
        inverse();
        expand(aligner, batch, /* ori */ 1);
        inverse();
    } else { /* ori = 0 */
        expand(aligner, batch, /* ori */ 1);
        expand(aligner, batch, /* ori */ -1);
    }
}

void Block::expand_end(PairAligner& aligner, int batch) {
    const std::vector<FragmentPtr>& fragments = fragments_;
    std::vector<int> main_end(fragments.size() - 1);
    std::vector<int> o_end(fragments.size() - 1);
    FragmentPtr main_f = fragments.back();
    while (true) {
        bool valid = true;
        BOOST_FOREACH (FragmentPtr f, fragments) {
            f->shift_end(batch);
            valid &= f->valid();
            f->shift_end(-batch);
        }
        if (!valid) {
            break;
        }
        std::string main_str = main_f->substr(-1, main_f->length() + batch);
        aligner.set_first(main_str.c_str(), main_str.size());
        for (int i = 0; i < fragments.size() - 1; i++) {
            FragmentPtr o_f = fragments[i];
            std::string o_str = o_f->substr(-1, o_f->length() + batch);
            aligner.set_second(o_str.c_str(), o_str.size());
            aligner.align(main_end[i], o_end[i]);
        }
        int min_end = *std::min_element(main_end.begin(), main_end.end());
        const float MIN_ACCEPTED = 0.5;
        if (min_end >= batch * MIN_ACCEPTED) {
            main_f->shift_end(min_end);
            for (int i = 0; i < fragments.size() - 1; i++) {
                FragmentPtr o_f = fragments[i];
                int delta = min_end - main_end[i];
                o_f->shift_end(o_end[i] - delta);
            }
        } else {
            break;
        }
    }
}

Block::Block()
{ }

std::ostream& operator<<(std::ostream& o, const Block& b) {
    int i = 0;
    BOOST_FOREACH (const FragmentPtr& f, b) {
        o << ">" << (++i) << std::endl; // FIXME
        o << *f << std::endl;
    }
    return o;
}

}

