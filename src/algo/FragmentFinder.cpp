/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "FragmentFinder.hpp"
#include "SeqI.hpp"
#include "Block.hpp"
#include "thread_pool.hpp"

namespace npge {

static bool pattern_length(const Processor* f,
                           std::string& message) {
    std::string p;
    p = f->opt_value("pattern").as<std::string>();
    Sequence::to_atgcn(p);
    if (p.empty()) {
        message = "'pattern' should not be empty";
        return false;
    }
    return true;
}

FragmentFinder::FragmentFinder() {
    add_opt("pattern", "Sequence searched for",
            std::string(""));
    add_opt_check(boost::bind(pattern_length, this, _1));
    add_gopt("max-matches", "Maximum number of matches",
             "MAX_MATCHES");
    add_opt_rule("max-matches >= 1");
    declare_bs("target", "search in sequences, add fragments");
}

class FinderTG : public ReusingThreadGroup,
    public SeqBase {
public:
    Fragments ff_;
    hash_t pattern_;
    std::string p_;
    size_t max_matches_;

    FinderTG(const FragmentFinder* f):
        SeqBase(*f->block_set()) {
        p_ = f->opt_value("pattern").as<std::string>();
        Sequence::to_atgcn(p_);
        anchor_ = p_.size();
        pattern_ = make_hash(p_.c_str(), p_.size(), 1);
        max_matches_ = f->opt_value("max-matches").as<int>();
        set_workers(f->workers());
        make_seqs();
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl();
};

class FinderWorker : public ThreadWorker {
public:
    Fragments ff_;

    FinderWorker(ThreadGroup* group):
        ThreadWorker(group) {
    }

    ~FinderWorker() {
        FinderTG* g = D_CAST<FinderTG*>(thread_group());
        Fragments& dst = g->ff_;
        dst.insert(dst.end(), ff_.begin(), ff_.end());
    }
};

class FinderTask : public ThreadTask, public SeqI {
public:
    Fragments& ff_;
    hash_t pattern_;
    std::string p_;
    size_t max_matches_;

    FinderTask(Sequence* seq, FinderWorker* w, FinderTG* g):
        ThreadTask(w),
        SeqI(seq, g),
        ff_(w->ff_),
        pattern_(g->pattern_),
        p_(g->p_),
        max_matches_(g->max_matches_) {
    }

    void add_match(int ori) {
        size_t min_pos = pos_;
        size_t max_pos = min_pos + anchor_ - 1;
        size_t begin = (ori == 1) ? min_pos : max_pos;
        if (anchor_ <= MAX_ANCHOR_SIZE ||
                seq_->substr(begin, anchor_, ori) == p_) {
            Fragment* f = new Fragment(seq_, min_pos,
                                       max_pos, ori);
            ff_.push_back(f);
        }
    }

    void test() {
        if (dir_ == pattern_) {
            add_match(1);
        }
        if (rev_ == pattern_) {
            add_match(-1);
        }
    }

    void run_impl() {
        init_state();
        test();
        size_t n = seq_->size() - anchor_;
        for (size_t i = 0; i < n; i++) {
            next_hash();
            test();
            if (ff_.size() > max_matches_) {
                return;
            }
        }
    }
};

ThreadTask* FinderTG::create_task_impl(ThreadWorker* worker) {
    if (it_ != end_) {
        Sequence* seq = *it_;
        it_++;
        FinderWorker* w = D_CAST<FinderWorker*>(worker);
        return new FinderTask(seq, w, this);
    } else {
        return 0;
    }
}

ThreadWorker* FinderTG::create_worker_impl() {
    return new FinderWorker(this);
}

void FragmentFinder::run_impl() const {
    FinderTG tg(this);
    tg.perform();
    BlockSet& bs = *block_set();
    int max_matches = opt_value("max-matches").as<int>();
    int matches = 0;
    BOOST_FOREACH (Fragment* f, tg.ff_) {
        if (matches < max_matches) {
            Block* block = new Block;
            block->insert(f);
            bs.insert(block);
        } else {
            delete f;
        }
    }
}

const char* FragmentFinder::name_impl() const {
    return "Locate fragment by its sequence";
}

}

