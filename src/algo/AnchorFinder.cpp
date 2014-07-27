/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <map>
#include "boost-xtime.hpp"
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include "po.hpp"

#include "AnchorFinder.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "BloomFilter.hpp"
#include "Exception.hpp"
#include "make_hash.hpp"
#include "complement.hpp"
#include "simple_task.hpp"
#include "thread_pool.hpp"
#include "throw_assert.hpp"
#include "SortedVector.hpp"
#include "boundaries.hpp"
#include "cast.hpp"

namespace npge {

AnchorFinder::AnchorFinder() {
    add_gopt("anchor-size", "anchor size", "ANCHOR_SIZE");
    add_gopt("anchor-fp",
             "Probability of false positive in Bloom filter "
             "(first step of AnchorFinder)", "ANCHOR_FP");
    add_opt_rule("anchor-size > 0");
    int max_anchor_size = sizeof(hash_t) * 8 / 2;
    add_opt_rule("anchor-size <= " + TO_S(max_anchor_size));
    declare_bs("target", "Blockset to search anchors in");
}

static int ns_in_fragment(const Fragment& f) {
    int result = 0;
    for (int i = 0; i < f.length(); i++) {
        if (f.raw_at(i) == 'N') {
            result += 1;
        }
    }
    return result;
}

struct CmpSeqSize {
    bool operator()(Sequence* a, Sequence* b) {
        return a->size() < b->size();
    }
};

typedef SortedVector<hash_t> Hashes;
typedef std::vector<Sequence*> Sequences;

struct AnchorFinderOptions {
    const AnchorFinder* finder_;
    BlockSet& bs_;

    typedef Sequences::iterator It;

    Sequences seqs_;
    It it_, end_;

    double error_prob_;
    int anchor_;

    AnchorFinderOptions(const AnchorFinder* f):
        finder_(f),
        bs_(*f->block_set()) {
        anchor_ = f->opt_value("anchor-size").as<int>();
        Decimal ep_d = f->opt_value("anchor-fp").as<Decimal>();
        error_prob_ = ep_d.to_d();
        //
        BOOST_FOREACH (const SequencePtr& s, bs_.seqs()) {
            seqs_.push_back(s.get());
        }
        // sort by size desc
        std::sort(seqs_.rbegin(), seqs_.rend(), CmpSeqSize());
        it_ = seqs_.begin();
        end_ = seqs_.end();
    }
};

class SeqI {
public:
    Sequence* seq_;
    size_t pos_;
    int ns_;
    int anchor_;

    hash_t dir_, rev_;

    SeqI(Sequence* seq, AnchorFinderOptions* opts):
        seq_(seq),
        pos_(0),
        ns_(0),
        anchor_(opts->anchor_) {
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
        // TODO why / 2
        return max_length / 2;
    } else {
        // TODO why * 1.5?
        return max_length / 2 * 3;
    }
}

class BloomTG : public ReusingThreadGroup,
    public AnchorFinderOptions {
public:
    BloomFilter bloom_;
    Hashes hashes_; // output

    BloomTG(const AnchorFinder* finder):
        AnchorFinderOptions(finder) {
        initialize_bloom();
    }

    static size_t pow4(int anchor_size) {
        ASSERT_LTE(anchor_size * 2, sizeof(size_t) * 8);
        return size_t(1) << (anchor_size * 2);
    }

    void initialize_bloom() {
        const BlockSet& bs = *(finder_->block_set());
        size_t length_sum = estimate_length(bs);
        if (anchor_ * 2 < sizeof(size_t) * 8) {
            size_t all_anchors = pow4(anchor_);
            if (all_anchors < length_sum) {
                length_sum = all_anchors;
            }
        }
        bloom_.set_members(length_sum, error_prob_);
        bloom_.set_optimal_hashes(length_sum);
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
    BloomFilter& bloom_;
    Hashes& hashes_;
    bool prev_;

    BloomTask(Sequence* seq, ThreadWorker* w):
        ThreadTask(w),
        SeqI(seq, D_CAST<BloomTG*>(thread_group())),
        bloom_(D_CAST<BloomTG*>(thread_group())->bloom_),
        hashes_(D_CAST<BloomWorker*>(worker())->hashes_) {
    }

    void test_and_add() {
        bool hash_found = false;
        if (ns_ == 0) {
            hash_t hash = std::min(dir_, rev_);
            hash_found = bloom_.test_and_add(hash);
            if (hash_found && !prev_) {
                hashes_.push_back(hash);
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

    FoundFragment(hash_t hash, Sequence* seq, size_t pos):
        hash_(hash), seq_(seq), pos_(pos) {
    }

    bool operator<(const FoundFragment& other) const {
        return hash_ < other.hash_;
    }

    bool operator>=(const FoundFragment& o) const {
        typedef boost::tuple<hash_t, Sequence*, size_t> Tie;
        return Tie(hash_, seq_, pos_) >=
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

static void fragmenttg_postprocess(FragmentTG& tg) {
    FFs& ffs = tg.ffs_;
    ffs.sort();
    ASSERT_TRUE(ffs.is_sorted_unique());
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
        }
    }
    if (block) {
        check_block(block, anchor);
    }
}

void AnchorFinder::run_impl() const {
    BloomTG bloomtg(this);
    bloomtg.perform();
    bloomtg_postprocess(bloomtg);
    FragmentTG fragmenttg(bloomtg.hashes_, this);
    fragmenttg.perform();
    bloomtg.hashes_.clear();
    fragmenttg_postprocess(fragmenttg);
}

const char* AnchorFinder::name_impl() const {
    return "Find anchors";
}

}

