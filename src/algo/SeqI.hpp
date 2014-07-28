/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#ifndef NPGE_SEQ_I_HPP_
#define NPGE_SEQ_I_HPP_

#include <algorithm>

#include "global.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "BlockSet.hpp"
#include "make_hash.hpp"
#include "complement.hpp"
#include "throw_assert.hpp"

namespace npge {

typedef std::vector<Sequence*> Sequences;

struct CmpSeqSize {
    bool operator()(Sequence* a, Sequence* b) {
        return a->size() < b->size();
    }
};

struct SeqBase {
    BlockSet& bs_;

    typedef Sequences::iterator It;

    Sequences seqs_;
    It it_, end_;

    int anchor_;

    SeqBase(BlockSet& bs):
        bs_(bs) {
    }

    void make_seqs() {
        seqs_.clear();
        BOOST_FOREACH (const SequencePtr& s, bs_.seqs()) {
            if (s->size() >= anchor_) {
                seqs_.push_back(s.get());
            }
        }
        // sort by size desc
        std::sort(seqs_.rbegin(), seqs_.rend(), CmpSeqSize());
        it_ = seqs_.begin();
        end_ = seqs_.end();
    }
};

inline int ns_in_fragment(const Fragment& f) {
    int result = 0;
    for (int i = 0; i < f.length(); i++) {
        if (f.raw_at(i) == 'N') {
            result += 1;
        }
    }
    return result;
}

class SeqI {
public:
    Sequence* seq_;
    size_t pos_;
    int ns_;
    int anchor_;

    hash_t dir_, rev_;

    SeqI(Sequence* seq, SeqBase* base):
        seq_(seq),
        pos_(0),
        ns_(0),
        anchor_(base->anchor_) {
    }

    void init_state() {
        ASSERT_GTE(seq_->size(), anchor_);
        Fragment init_f(seq_, 0, anchor_ - 1);
        ns_ = ns_in_fragment(init_f);
        dir_ = init_f.hash();
        init_f.inverse();
        rev_ = init_f.hash();
        ASSERT_EQ(rev_, complement_hash(dir_, anchor_));
    }

    void update_hash(hash_t& hash, char remove_char,
                     char add_char, bool direct) {
        hash = reuse_hash(hash, anchor_,
                          remove_char, add_char,
                          direct);
    }

    void next_hash() {
        char remove_char = seq_->char_at(pos_);
        char add_char = seq_->char_at(pos_ + anchor_);
        pos_ += 1;
        if (remove_char == 'N') {
            ns_ -= 1;
        }
        if (add_char == 'N') {
            ns_ += 1;
        }
        update_hash(dir_, remove_char, add_char, true);
        remove_char = complement(remove_char);
        add_char = complement(add_char);
        update_hash(rev_, remove_char, add_char, false);
    }
};

}

#endif

