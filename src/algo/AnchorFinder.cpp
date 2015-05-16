/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2012-2015 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include "boost-xtime.hpp"
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "AnchorFinder.hpp"
#include "SeqI.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "BloomFilter.hpp"
#include "Exception.hpp"
#include "thread_pool.hpp"
#include "throw_assert.hpp"
#include "SortedVector.hpp"
#include "boundaries.hpp"
#include "cast.hpp"

namespace npge {

typedef SortedVector<hash_t> Hashes;

struct AnchorFinderImpl {
    Hashes used_hashes_;
};

struct AnchorFinder::Impl : public AnchorFinderImpl {
};

AnchorFinder::AnchorFinder():
    impl_(new Impl) {
    add_gopt("anchor-size", "anchor size", "ANCHOR_SIZE");
    add_gopt("anchor-fp",
             "Probability of false positive in Bloom filter "
             "(first step of AnchorFinder)", "ANCHOR_FP");
    add_opt("anchor-similar",
            "If neighbour anchors are skipped",
            true);
    add_gopt("max-anchor-fragments",
             "Maximum number of anchors fragments to return",
             "MAX_ANCHOR_FRAGMENTS");
    add_opt_rule("anchor-size > 0");
    int max_anchor_size = sizeof(hash_t) * 8 / 2;
    add_opt_rule("anchor-size <= " + TO_S(MAX_ANCHOR_SIZE));
    declare_bs("target", "Blockset to search anchors in");
}

AnchorFinder::~AnchorFinder() {
    delete impl_;
}

struct AnchorFinderOptions : public SeqBase {
    const AnchorFinder* finder_;

    double error_prob_;
    int max_anchor_fragments_;
    bool similar_;

    AnchorFinderOptions(const AnchorFinder* f):
        SeqBase(*f->block_set()),
        finder_(f) {
        anchor_ = f->opt_value("anchor-size").as<int>();
        Decimal ep_d = f->opt_value("anchor-fp").as<Decimal>();
        error_prob_ = ep_d.to_d();
        similar_ = f->opt_value("anchor-similar").as<bool>();
        max_anchor_fragments_ =
            f->opt_value("max-anchor-fragments").as<int>();
        make_seqs();
    }
};

// Bloom filter

static size_t estimate_length(const BlockSet& bs) {
    typedef std::map<std::string, size_t> GenomeToLength;
    GenomeToLength gtl;
    size_t max_length = 0;
    BOOST_FOREACH (SequencePtr s, bs.seqs()) {
        size_t& genome_length = gtl[s->genome()];
        genome_length += s->size();
        max_length = std::max(max_length, genome_length);
    }
    if (gtl.size() == 1) {
        // consensuses
        return max_length;
    } else {
        // TODO why * 1.5?
        return max_length / 2 * 3;
    }
}

class BloomTG : public ReusingThreadGroup,
    public AnchorFinderOptions {
public:
    const Hashes& used_;
    BloomFilter bloom_;
    Hashes hashes_; // output
    size_t length_sum_;

    BloomTG(const AnchorFinder* finder,
            const Hashes& used_hashes):
        AnchorFinderOptions(finder),
        used_(used_hashes) {
        set_workers(finder->workers());
        initialize_bloom();
    }

    static size_t pow4(int anchor_size) {
        ASSERT_LTE(anchor_size * 2, sizeof(size_t) * 8);
        return size_t(1) << (anchor_size * 2);
    }

    void initialize_bloom() {
        const BlockSet& bs = *(finder_->block_set());
        length_sum_ = estimate_length(bs);
        if (anchor_ * 2 < sizeof(size_t) * 8) {
            size_t all_anchors = pow4(anchor_);
            if (all_anchors < length_sum_) {
                length_sum_ = all_anchors;
            }
        }
        bloom_.set_members(length_sum_, error_prob_);
        bloom_.set_optimal_hashes(length_sum_);
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl();
};

class BloomWorker : public ThreadWorker {
public:
    Hashes hashes_;

    BloomWorker(ThreadGroup* group):
        ThreadWorker(group) {
    }

    ~BloomWorker() {
        BloomTG* g = D_CAST<BloomTG*>(thread_group());
        g->hashes_.extend(hashes_);
    }
};

class BloomTask : public ThreadTask, public SeqI {
public:
    const Hashes& used_;
    BloomFilter& bloom_;
    Hashes& hashes_;
    bool prev_;
    bool similar_;

    BloomTask(Sequence* seq, ThreadWorker* w):
        ThreadTask(w),
        SeqI(seq, D_CAST<BloomTG*>(thread_group())),
        used_(D_CAST<BloomTG*>(thread_group())->used_),
        bloom_(D_CAST<BloomTG*>(thread_group())->bloom_),
        hashes_(D_CAST<BloomWorker*>(worker())->hashes_),
        prev_(false),
        similar_(D_CAST<BloomTG*>(thread_group())->similar_) {
    }

    void test_and_add() {
        bool hash_found = false;
        if (ns_ == 0) {
            hash_t hash = std::min(dir_, rev_);
            if (!used_.has_elem(hash)) {
                hash_found = bloom_.test_and_add(hash);
                if (hash_found && (!prev_ || !similar_)) {
                    hashes_.push_back(hash);
                }
            }
        }
        prev_ = hash_found;
    }

    void run_impl() {
        if (seq_->size() < anchor_) {
            return;
        }
        init_state();
        prev_ = false;
        test_and_add();
        size_t n = seq_->size() - anchor_;
        for (size_t i = 0; i < n; i++) {
            next_hash();
            test_and_add();
        }
    }
};

ThreadTask* BloomTG::create_task_impl(ThreadWorker* worker) {
    if (it_ != end_) {
        Sequence* seq = *it_;
        it_++;
        return new BloomTask(seq, worker);
    } else {
        return 0;
    }
}

ThreadWorker* BloomTG::create_worker_impl() {
    return new BloomWorker(this);
}

static void bloomtg_postprocess(BloomTG& g) {
    g.bloom_.clear();
    Hashes& hashes = g.hashes_;
    hashes.sort();
    hashes.unique();
}

// find fragments matching found hashes

struct FoundFragment {
    hash_t hash_;
    Sequence* seq_;
    size_t pos_; // if ori = -1, pos = size + min_pos

    FoundFragment() {
    }

    FoundFragment(hash_t hash, Sequence* seq, size_t pos):
        hash_(hash), seq_(seq), pos_(pos) {
    }

    bool operator<(const FoundFragment& o) const {
        typedef boost::tuple<hash_t, Sequence*, size_t> Tie;
        return Tie(hash_, seq_, pos_) <
               Tie(o.hash_, o.seq_, o.pos_);
    }

    Fragment* make_fragment(int anchor) const {
        bool direct = (pos_ < seq_->size());
        int ori = direct ? 1 : -1;
        size_t min_pos = direct ? pos_ : (pos_ - seq_->size());
        size_t max_pos = min_pos + anchor - 1;
        return new Fragment(seq_, min_pos, max_pos, ori);
    }
};

typedef SortedVector<FoundFragment> FFs;

class FragmentTG : public ReusingThreadGroup,
    public AnchorFinderOptions {
public:
    const Hashes& hashes_; // input
    FFs ffs_; // output

    FragmentTG(const Hashes& hashes,
               const AnchorFinder* finder):
        AnchorFinderOptions(finder),
        hashes_(hashes) {
        set_workers(finder->workers());
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl();
};

class FragmentWorker : public ThreadWorker {
public:
    FFs ffs_;

    FragmentWorker(ThreadGroup* group):
        ThreadWorker(group) {
    }

    ~FragmentWorker() {
        FragmentTG* g = D_CAST<FragmentTG*>(thread_group());
        g->ffs_.extend(ffs_);
    }
};

class FragmentTask : public ThreadTask, public SeqI {
public:
    const Hashes& hashes_; // input
    FFs& ffs_; // output

    FragmentTask(Sequence* seq, ThreadWorker* w):
        ThreadTask(w),
        SeqI(seq, D_CAST<FragmentTG*>(thread_group())),
        hashes_(D_CAST<FragmentTG*>(thread_group())->hashes_),
        ffs_(D_CAST<FragmentWorker*>(worker())->ffs_) {
    }

    void push(hash_t hash, bool direct) {
        size_t pos = pos_;
        if (direct == false) {
            pos += seq_->size();
        }
        ffs_.push_back(FoundFragment(hash, seq_, pos));
    }

    void test_and_push() {
        if (ns_ == 0) {
            hash_t hash = std::min(dir_, rev_);
            bool hash_found = hashes_.has_elem(hash);
            if (hash_found) {
                bool direct = (hash == dir_);
                push(hash, direct);
            }
        }
    }

    void run_impl() {
        if (seq_->size() < anchor_) {
            return;
        }
        init_state();
        test_and_push();
        size_t n = seq_->size() - anchor_;
        for (size_t i = 0; i < n; i++) {
            next_hash();
            test_and_push();
        }
    }
};

ThreadTask* FragmentTG::create_task_impl(ThreadWorker* worker) {
    if (it_ != end_) {
        Sequence* seq = *it_;
        it_++;
        return new FragmentTask(seq, worker);
    } else {
        return 0;
    }
}

ThreadWorker* FragmentTG::create_worker_impl() {
    return new FragmentWorker(this);
}

static void check_block(const Block* block, int anchor) {
    int size = block->size();
    ASSERT_GTE(size, 2);
    int length = block->alignment_length();
    ASSERT_EQ(length, anchor);
    Fragments ff(block->begin(), block->end());
    for (int pos = 0; pos < length; pos++) {
        char c = ff[0]->raw_at(pos);
        for (int f_i = 1; f_i < size; f_i++) {
            ASSERT_EQ(ff[f_i]->raw_at(pos), c);
        }
    }
}

static void fragmenttg_postprocess(FragmentTG& tg,
                                   Hashes& used_hashes) {
    FFs& ffs = tg.ffs_;
    ffs.sort();
    ASSERT_TRUE(ffs.is_sorted_unique());
    if (ffs.size() > tg.max_anchor_fragments_) {
        ffs.resize(tg.max_anchor_fragments_);
    }
    int anchor = tg.anchor_;
    BlockSet& bs = tg.bs_;
    const FoundFragment* prev = 0;
    Block* block = 0;
    BOOST_FOREACH (const FoundFragment& ff, ffs) {
        if (prev == 0) {
            prev = &ff;
        } else if (prev->hash_ == ff.hash_) {
            if (block == 0) {
                block = new Block;
                bs.insert(block);
                block->insert(prev->make_fragment(anchor));
            }
            ASSERT_TRUE(block);
            block->insert(ff.make_fragment(anchor));
        } else {
            prev = &ff;
            if (block) {
                check_block(block, anchor);
            }
            block = 0;
            used_hashes.push_back(ff.hash_);
        }
    }
    if (block) {
        check_block(block, anchor);
    }
}

void AnchorFinder::run_impl() const {
    BloomTG bloomtg(this, impl_->used_hashes_);
    bloomtg.perform();
    bloomtg_postprocess(bloomtg);
    FragmentTG fragmenttg(bloomtg.hashes_, this);
    fragmenttg.perform();
    bloomtg.hashes_.clear();
    bool sort_used_hashes = !impl_->used_hashes_.empty();
    fragmenttg_postprocess(fragmenttg, impl_->used_hashes_);
    if (sort_used_hashes) {
        impl_->used_hashes_.sort();
    }
    ASSERT_TRUE(impl_->used_hashes_.is_sorted_unique());
}

const char* AnchorFinder::name_impl() const {
    return "Find anchors";
}

}

