/*
 * bloomrepeats, Find genomic repeats, using Bloom filter based prefiltration
 * Copyright (C) 2012 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

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

}

