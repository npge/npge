/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2013 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <climits>
#include <vector>
#include <set>
#include <string>
#include <boost/cast.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>

#include "block_hash.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"
#include "AlignmentRow.hpp"
#include "Block.hpp"
#include "BlockSet.hpp"
#include "thread_pool.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"
#include "global.hpp"

namespace npge {

hash_t block_hash(const Block* block) {
    Strings ids_dir, ids_inv;
    BOOST_FOREACH (Fragment* f, *block) {
        ids_dir.push_back(f->id());
        f->inverse(/* inverse_row */ false);
        ids_inv.push_back(f->id());
        f->inverse(/* inverse_row */ false);
    }
    std::sort(ids_dir.begin(), ids_dir.end());
    std::sort(ids_inv.begin(), ids_inv.end());
    Strings& ids = (ids_dir < ids_inv) ? ids_dir : ids_inv;
    std::string joint = boost::algorithm::join(ids, " ");
    const int LOOP_SIZE = sizeof(hash_t) * 2;
    // 2 = for * and for ^
    int new_size = ((joint.size() + LOOP_SIZE - 1) / LOOP_SIZE)
                   * LOOP_SIZE;
    joint.resize(new_size, ' ');
    const char* c = joint.c_str();
    const hash_t* value = reinterpret_cast<const hash_t*>(c);
    int loops = joint.size() / LOOP_SIZE;
    hash_t a = 1;
    for (int i = 0; i < loops; i++) {
        a *= value[2 * i];
        a ^= value[2 * i + 1];
    }
    return a;
}

class HashTask;
class HashWorker;
class HashGroup;

class HashGroup : public ReusingThreadGroup {
public:
    HashGroup(const BlockSet& block_set):
        it_(block_set.begin()),
        end_(block_set.end()), hash_(0) {
    }

    ThreadTask* create_task_impl(ThreadWorker* worker);

    ThreadWorker* create_worker_impl();

    BlockSet::const_iterator it_;
    BlockSet::const_iterator end_;

    hash_t hash_;
};

class HashWorker : public ThreadWorker {
public:
    HashWorker(HashGroup* thread_group):
        ThreadWorker(thread_group), hash_(0) {
    }

    ~HashWorker() {
        HashGroup* g = D_CAST<HashGroup*>(thread_group());
        g->hash_ ^= hash_;
    }

    hash_t hash_;
};

class HashTask : public ThreadTask {
public:
    HashTask(const Block* block, HashWorker* worker):
        ThreadTask(worker), block_(block) {
    }

    void run_impl() {
        HashWorker* w = D_CAST<HashWorker*>(worker());
        w->hash_ ^= block_hash(block_);
    }

private:
    const Block* block_;
};

ThreadTask* HashGroup::create_task_impl(ThreadWorker* worker) {
    while (it_ != end_ && (*it_)->size() <= 1) {
        it_++;
    }
    if (it_ == end_) {
        return 0;
    } else {
        HashWorker* w = D_CAST<HashWorker*>(worker);
        HashTask* task = new HashTask(*it_, w);
        it_++;
        return task;
    }
}

ThreadWorker* HashGroup::create_worker_impl() {
    return new HashWorker(this);
}

hash_t blockset_hash(const BlockSet& block_set, int workers) {
    HashGroup hash_group((block_set));
    hash_group.set_workers(workers);
    hash_group.perform();
    return hash_group.hash_;
}

std::string block_id(const Block* block) {
    return TO_S(block->size()) + "x" +
           TO_S(block->alignment_length());
}

typedef std::set<std::string> StringSet;

bool has_repeats(const Block* block) {
    StringSet genomes;
    BOOST_FOREACH (Fragment* f, *block) {
        if (f->seq()) {
            std::string genome = f->seq()->genome();
            if (genomes.find(genome) != genomes.end()) {
                return true;
            }
            genomes.insert(genome);
        }
    }
    return false;
}

bool is_exact_stem(const Block* block, int genomes) {
    return block->size() == genomes && !has_repeats(block);
}

int genomes_number(const BlockSet& block_set) {
    StringSet all_genomes;
    BOOST_FOREACH (const SequencePtr& seq, block_set.seqs()) {
        all_genomes.insert(seq->genome());
    }
    return all_genomes.size();
}

std::string block_name(const Block* b, int genomes) {
    char type = (b->size() == 1) ? 'u' :
                has_repeats(b) ? 'r' :
                (b->size() == genomes) ? 's' : 'h';
    std::string name;
    name += type;
    name += block_id(b);
    return name;
}

bool has_alignment(const Block* block) {
    BOOST_FOREACH (const Fragment* f, *block) {
        if (!f->row()) {
            return false;
        }
    }
    return true;
}

static Fragment* min_fragment(const Block* block) {
    ASSERT_GTE(block->size(), 1);
    Fragment* result = block->front();
    BOOST_FOREACH (Fragment* f, *block) {
        if (*f < *result) {
            result = f;
        }
    }
    return result;
}

bool block_less(const Block* a, const Block* b) {
    int a_size = a->size();
    int b_size = b->size();
    if (a_size < b_size) {
        return true;
    } else if (a_size > b_size) {
        return false;
    } else if (a_size == 0) {
        return false;
    }
    int a_length = a->alignment_length();
    int b_length = b->alignment_length();
    if (a_length < b_length) {
        return true;
    } else if (a_length > b_length) {
        return false;
    }
    Fragment* a_f = min_fragment(a);
    Fragment* b_f = min_fragment(b);
    return *a_f < *b_f;
}

bool block_greater(const Block* a, const Block* b) {
    if (block_less(b, a)) {
        return true;
    } else {
        return false;
    }
}

static void test_fragment(Fragment* f, int length) {
    ASSERT_LTE(f->length(), f->row()->length());
    ASSERT_MSG(f->row()->length() == length,
               ("Length of row of fragment " + f->id() +
                " (" + TO_S(f->row()->length()) + ") "
                "differs from block alignment length"
                " (" + TO_S(length) + ")").c_str());
}

void test_block(const Block* block) {
    if (block->empty()) {
        return;
    }
    bool has_rows = (block->front()->row() != 0);
    int length = block->alignment_length();
    BOOST_FOREACH (Fragment* f, *block) {
        bool f_has_rows = (f->row() != 0);
        ASSERT_EQ(f_has_rows, has_rows);
        if (f_has_rows) {
            test_fragment(f, length);
        }
    }
}

}

